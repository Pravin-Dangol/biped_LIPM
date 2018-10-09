% simulate the leg mechanism
% Notes:
%   The leg-end is steered to several way points as the IK problem is resolved 
%   using 1) a nonlinear solver "fsolve", and 2) a trigonometric-based analytical
%   method.
%
%   q1: HD1 angle
%   q2: HD2 angle
%   alpha1: (other leg joint angles, refer to my notes)
%   alpha2:
%   alpha3:
%   alpha4:
%   alpha5:
%   alpha6:
%
% 12/22/2017 Writen by Alirza Ramezani, Pasadena, CA
clear all; close all; clc;


% param: from CAD model
% R1 R2 R3 Rc1 Rc2 Rc3 e1 e2 e3 e4 [mm]
func_params;

% desired leg-end trajectory
stp_ln = -0.05;
stp_ht = 0.1;
n = 20; % number of way points
p_traj = func_bezier(stp_ln,stp_ht,n);
diff_p_traj = diff(p_traj);
diff_p_traj = [diff_p_traj;[0,0]];

q_traj = []; % stack HD angles (e.g., [q1,q2])
alpha_sol_num = [];    % stack unknown angles (e.g., alpha1, alpha2, etc.)
alpha_sol_analytical = [];    % stack unknown angles (e.g., alpha1, alpha2, etc.)
joint_pos_num = [];   % stack joint positions...
joint_pos_analytical = [];   % stack joint positions...

% these are the exact solutions (alphas) saved at "func_params" from the
% solver for LEO's leg when q1 = 99.99 [deg] and q2 =159.99 [deg] 
xi = x0;

% HD angles
qi = [q1,q2].';

for i=1:n
    
    disp(['determining mechanism config. for number ', num2str(i),' way-point out of ',num2str(n),' way-points ... please wait...']);
    
    % solve IK problem
    
    % commanded movements 
    delta_p = diff_p_traj(i,:).';  % for not just move straight ...
    
    % input:
    %   q(1): HD1 angle
    %   q(2): HD2 angle
    %
    %   qp(1): alpha1 (look at notes)
    %   qp(2): alpha2
    %   qp(3): alpha3
    %   qp(4): alpha4
    %   qp(5): alpha5
    %   qp(6): alpha6
    %
    %   delta_p(1): horizontal variations in leg-end position
    %   delta_p(2): vertical variations in leg-end position
    %
    %   params: R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4
    %
    % output:
    %   delta_q(1): variations in HD1 angle
    %   delta_q(2): variations in HD1 angle
    %
    [delta_q] = func_task_space_to_joint_space(qi,xi,delta_p,R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4);
    
    % update desired HD angles
    qi = qi + delta_q;
    
    % solve for unknowns....
    % x_sol: [rad]
    % alpha1
    % alpha2
    % alpha3
    % alpha4
    % alpha5
    % alpha6
    
    % numerical solution
    [x_sol_num]=func_sol_mechanism(xi,qi(1),qi(2),R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4);
    xi = x_sol_num.'; % Note: I need to update the initial guess xi for the the numeric solver
    
    % analytical solution
    [a1,a2,a3,a4,a5,a6]=func_q_prime_trigonometry(qi(1), qi(2), R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4);
    x_sol_analytical=[a1,a2,a3,a4,a5,a6].';
    
    % computing joint positions ...
    % Note: these are symbolically generated...
    
    disp('computing fwd kinematics terms ...');
    % computing fwd kinematics terms ...
    % Note: these are symbolically generated...
    % x:
    % alpha1
    % alpha2
    % alpha3
    % alpha4
    % alpha5
    % alpha6
    % motor angles:
    % q1
    % q2
    % params:
    % R1 R2 R3 Rc1 Rc2 Rc3 e1 e2 e3 e4
    
    % output:
    % p1, p2 ,p3, ... p9
    % v1, v2 ,v3, ... v9
    [p1,p2,p3,p4,p5,p6,p7,p8,p9,v1,v2,v3,v4,v5,v6,v7,v8,v9] = ...
        func_leg_fwd_kin(x_sol_num, qi(1), qi(2), R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4);
    joint_pos_num = [joint_pos_num;p1.',p2.',p3.',p4.',p5.',p6.',p7.',p8.',p9.'];
    
    [p1,p2,p3,p4,p5,p6,p7,p8,p9,v1,v2,v3,v4,v5,v6,v7,v8,v9] = ...
    func_leg_fwd_kin(x_sol_analytical, qi(1), qi(2), R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4);
    joint_pos_analytical = [joint_pos_analytical;p1.',p2.',p3.',p4.',p5.',p6.',p7.',p8.',p9.'];
    
    
    % stack data
    q_traj = [q_traj;qi.'];
    alpha_sol_num = [alpha_sol_num;x_sol_num.'];
    alpha_sol_analytical = [alpha_sol_analytical;x_sol_analytical.'];
    
end

q_traj = q_traj*180/pi;
alpha_sol_num = alpha_sol_num*180/pi;
alpha_sol_analytical = alpha_sol_analytical*180/pi;

% plot results ...
%%
scrsz = get(groot,'ScreenSize');
fh = figure('Name','mechanism angles',...
    'Renderer','opengl',...
    'GraphicsSmoothing','on');
ah = axes('Box','on',...
    'XGrid','off',...
    'YGrid','off',...
    'Parent',fh);
xlabel(ah,'[time steps]');
ylabel(ah,'[deg]');
hold(ah,'on');
plot(ah,1:n,q_traj(:,1),'-k');
plot(ah,1:n,q_traj(:,2),'--k');
plot(ah,1:n,alpha_sol_num(:,1),'-b');
plot(ah,1:n,alpha_sol_num(:,2),'--b');
plot(ah,1:n,alpha_sol_num(:,3),'-r');
plot(ah,1:n,alpha_sol_num(:,4),'--r');
plot(ah,1:n,alpha_sol_num(:,5),'-g');
plot(ah,1:n,alpha_sol_num(:,6),'--g');
plot(ah,1:n,alpha_sol_analytical(:,1),'-b','LineWidth',2);
plot(ah,1:n,alpha_sol_analytical(:,2),'--b','LineWidth',2);
plot(ah,1:n,alpha_sol_analytical(:,3),'-r','LineWidth',2);
plot(ah,1:n,alpha_sol_analytical(:,4),'--r','LineWidth',2);
plot(ah,1:n,alpha_sol_analytical(:,5),'-g','LineWidth',2);
plot(ah,1:n,alpha_sol_analytical(:,6),'--g','LineWidth',2);
legend(ah,'q_1','q_2','\alpha_1^{num}','\alpha_2^{num}','\alpha_3^{num}','\alpha_4^{num}','\alpha_5^{num}','\alpha_6^{num}'...
    ,'\alpha_1^{an}','\alpha_2^{an}','\alpha_3^{an}','\alpha_4^{an}','\alpha_5^{an}','\alpha_6^{an}');

%%
scrsz = get(groot,'ScreenSize');
fh = figure('Name','Leg-end trajectory: analytical versus numerical solution',...
    'Renderer','opengl',...
    'GraphicsSmoothing','on');
ah = axes('Box','on',...
    'XGrid','off',...
    'YGrid','off',...
    'DataAspectRatio',[1,1,1],...
    'PlotBoxAspectRatio',[1,1,1],...
    'Parent',fh);
xlabel(ah,'[m]');
ylabel(ah,'[m]');
hold(ah,'on');

plot(ah,joint_pos_num(:,17),joint_pos_num(:,18),'o-k');
plot(ah,joint_pos_analytical(:,17),joint_pos_analytical(:,18),'o-r');
legend(ah,'numeric','analytic');



%%
% animate mechanism ...
pos = joint_pos_num;

scrsz = get(groot,'ScreenSize');
fh = figure('Name','LEO leg mechanism (side view)',...
    'Renderer','opengl',...
    'GraphicsSmoothing','on');
ah = axes('Box','on',...
    'XGrid','off',...
    'YGrid','off',...
    'DataAspectRatio',[1,1,1],...
    'PlotBoxAspectRatio',[1,1,1],...
    'Parent',fh);
xlabel(ah,'[m]');
ylabel(ah,'[m]');
ylim(ah,[-0.5 0]);
hold(ah,'off');


for i=1:n
 
    cla(ah);    
 
    p1 = pos(i,1:2);
    p2 = pos(i,3:4);
    p3 = pos(i,5:6);
    p4 = pos(i,7:8);
    p5 = pos(i,9:10);
    p6 = pos(i,11:12);
    p7 = pos(i,13:14);
    p8 = pos(i,15:16);
    p9 = pos(i,17:18);
    
    % coordinates of the virtual leg
    p_h = (p1+p2)./2;
    p_e = p9;
    
    % plot links ...
    
    % first loop
    line(ah,[p1(1),p2(1)],[p1(2),p2(2)],'Color','b','Marker','o');
    line(ah,[p2(1),p3(1)],[p2(2),p3(2)],'Color','b','Marker','o');
    line(ah,[p3(1),p4(1)],[p3(2),p4(2)],'Color','b','Marker','o');
    line(ah,[p4(1),p5(1)],[p4(2),p5(2)],'Color','b','Marker','o');
    line(ah,[p5(1),p1(1)],[p5(2),p1(2)],'Color','b','Marker','o');
    line(ah,[p1(1),p3(1)],[p1(2),p3(2)],'Color','r','LineStyle','--');  
    line(ah,[p3(1),p5(1)],[p3(2),p5(2)],'Color','r','LineStyle','--');  
    
    % second loop
    line(ah,[p6(1),p4(1)],[p6(2),p4(2)],'Color','b','Marker','o');
    line(ah,[p7(1),p6(1)],[p7(2),p6(2)],'Color','b','Marker','o');
    line(ah,[p8(1),p7(1)],[p8(2),p7(2)],'Color','b','Marker','o');
    line(ah,[p4(1),p8(1)],[p4(2),p8(2)],'Color','b','Marker','o');
    line(ah,[p9(1),p7(1)],[p9(2),p7(2)],'Color','b','Marker','o');
    line(ah,[p6(1),p8(1)],[p6(2),p8(2)],'Color','r','LineStyle','--');  

    % virtual leg
    line(ah,[p_h(1),p_e(1)],[p_h(2),p_e(2)],'Color','g','LineStyle','--');  
    
    % 
    ah.DataAspectRatio = [1,1,1];
    ah.PlotBoxAspectRatio = [1,1,1];
        
    %
    drawnow;
    pause(0.02);
end

%%



































