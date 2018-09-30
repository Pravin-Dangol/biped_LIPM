% simulate the leg mechanism
% Notes:
%   1) q1 (HD1) and q2 (HD2) are given and a nonlinear solver "fsolve" is used to 
%   identify the configuration of LEO's leg, which is made out of two
%   loops. The leg has 2-DoFs. The angles are:
%   
%   q1:
%   q2:
%   alpha1:
%   alpha2:
%   alpha3:
%   alpha4:
%   alpha5:
%   alpha6:
%
% 12/22/2017 Writen by Alirza Ramezani, Pasadena, CA
clear all; close all; clc;


%param: from CAD model
% R1 R2 R3 Rc1 Rc2 Rc3 e1 e2 e3 e4 [mm]
func_params;

% Note: if swing_mode flag is set:
%       true - in a ramp increases the q1 & q2 angles to generate swing motion...
%       false - in a ramp increases q1  and decreases q2 angles to collaps
%       the leg mechanism...
swing_mode = true;

delta = 20; % change in q1 & q2 angles at the end of the ramp [deg]
step_sz = 1; % change in q1 & q2 angle at each time step [deg]

if swing_mode
    q1_tmp = q1:step_sz*pi/180:q1+delta*pi/180;     % q1 positions...
    q2_tmp = q2:step_sz*pi/180:q2+delta*pi/180;     % q2 positions...
else
    q1_tmp = q1:step_sz*pi/180:q1+delta*pi/180;     % q1 positions...
    q2_tmp = q2:-step_sz*pi/180:q2-delta*pi/180;     % q2 positions...
end

input = [q1_tmp.',q2_tmp.'];

output = [];    % stack unknown angles (e.g., alpha1, alpha2, etc.)
pos = [];   % stack joint positions...


n = length(input);

xi = x0;

for i=1:n
    
    disp(['determining mechanism config. for number ', num2str(i),' way-point out of ',num2str(n),' way-points ... please wait...']);
    
    q1 = input(i,1);
    q2 = input(i,2);
    
    % solve for unknowns....
    % x_sol: [rad]
    % alpha1
    % alpha2
    % alpha3
    % alpha4
    % alpha5
    % alpha6
    [x_sol]=func_sol_mechanism(xi,q1,q2,R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4);
    
    xi = x_sol.';
    
    alpha1 = x_sol(1);
    alpha2 = x_sol(2);
    alpha3 = x_sol(3);
    alpha4 = x_sol(4);
    alpha5 = x_sol(5);
    alpha6 = x_sol(6);
    
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
    func_leg_fwd_kin(x_sol, q1, q2, R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4);
    
    % stacking in row...
    pos = [pos;p1.',p2.',p3.',p4.',p5.',p6.',p7.',p8.',p9.'];
    output = [output;x_sol.'];
    
end

input = input*180/pi;
output = output*180/pi;

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
plot(ah,1:n,input(:,1),'-k');
plot(ah,1:n,input(:,2),'--k');
plot(ah,1:n,output(:,1),'-b');
plot(ah,1:n,output(:,2),'--b');
plot(ah,1:n,output(:,3),'-r');
plot(ah,1:n,output(:,4),'--r');
plot(ah,1:n,output(:,5),'-g');
plot(ah,1:n,output(:,6),'--g');
legend(ah,'q_1','q_2','\alpha_1','\alpha_2','\alpha_3','\alpha_4','\alpha_5','\alpha_6');

%%
% animate mechanism ...
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
hold(ah,'on');

for i=1:n
    p1 = pos(i,1:2); 
    p2 = pos(i,3:4);
    p3 = pos(i,5:6);
    p4 = pos(i,7:8);
    p5 = pos(i,9:10);
    p6 = pos(i,11:12);
    p7 = pos(i,13:14);
    p8 = pos(i,15:16);
    p9 = pos(i,17:18);
    
    % plot joints ...
    plot(ah,p1(1),p1(2),'ko');
    plot(ah,p2(1),p2(2),'ko');
    plot(ah,p3(1),p3(2),'ko');
    plot(ah,p4(1),p4(2),'ko');
    plot(ah,p5(1),p5(2),'ko');
    plot(ah,p6(1),p6(2),'ko');
    plot(ah,p7(1),p7(2),'ko');
    plot(ah,p8(1),p8(2),'ko');
    plot(ah,p9(1),p9(2),'ko');
    
    % plot links ...
    % first loop
    line(ah,[p1(1),p2(1)],[p1(2),p2(2)],'Color','b');
    line(ah,[p2(1),p3(1)],[p2(2),p3(2)],'Color','b');
    line(ah,[p3(1),p4(1)],[p3(2),p4(2)],'Color','b');
    line(ah,[p4(1),p5(1)],[p4(2),p5(2)],'Color','b');
    line(ah,[p5(1),p1(1)],[p5(2),p1(2)],'Color','b');
    
    % second loop
    line(ah,[p6(1),p4(1)],[p6(2),p4(2)],'Color','b');
    line(ah,[p7(1),p6(1)],[p7(2),p6(2)],'Color','b');
    line(ah,[p8(1),p7(1)],[p8(2),p7(2)],'Color','b');
    line(ah,[p4(1),p8(1)],[p4(2),p8(2)],'Color','b');
    line(ah,[p9(1),p7(1)],[p9(2),p7(2)],'Color','b');
    
    drawnow;
    pause(1);
end





































