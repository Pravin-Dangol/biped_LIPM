% optimize leg dims defined in "func_params.m" subject to a cost function:
%   q1
%   q2
%   R1
%   R2
%   R3
%   Rc1
%   Rc2
%   Rc3
%   e1
%   e2
%   e3
%   e4
%   alpha1
%   alpha2
%   alpha3
%   alpha4
%   alpha5
%   alpha6
%   w: LEO's weight
%
%   x_sol:
%       q1
%       q2
%       R1
%       R2
%       R3
%       Rc1
%       Rc2
%       Rc3
%       e1
%       e2
%       e3
%       e4
%   fval,exitflag,output,lambda,grad,hessian
% by Alireza Ramezani, 5-19-2018, Pasadena, CA
function [x_sol,fval,exitflag,output,lambda,grad,hessian] = func_opt_leg_dims(q1, q2, R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, w)


% set opt settings...
% options = optimoptions('fmincon',...
%     'Display','iter',...
%     'Algorithm','sqp',...
%     'ConstraintTolerance',1e-4,...
%     'FunctionTolerance', 1,...
%     'HonorBounds',true,...
%     'MaxFunctionEvaluations',1e3,...
%     'MaxIterations',1e3,...
%     'ScaleProblem',false,...
%     'UseParallel',false);

options = optimoptions('fmincon',...
    'Display','iter',...
    'Algorithm','sqp',...
    'ConstraintTolerance',1e-4,...
    'FunctionTolerance', 1,...
    'HonorBounds',true,...
    'MaxFunctionEvaluations',1e3,...
    'MaxIterations',1e3,...
    'ScaleProblem','none',...
    'UseParallel',false);

%params
% func_params;

% Note: ????
alpha0 = zeros(1,6);
alpha0(1) = alpha1;
alpha0(2) = alpha2;
alpha0(3) = alpha3;
alpha0(4) = alpha4;
alpha0(5) = alpha5;
alpha0(6) = alpha6;

%   x0:
%       q1
%       q2
%       R1
%       R2
%       R3
%       Rc1
%       Rc2
%       Rc3
%       e1
%       e2
%       e3
%       e4
x0 = [q1, q2, R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4].';

% global vars
GRF = [0, w].'; % ground reaction force defined in "func_params.m"
alpha_sol = [];

% joint pos. and vel.
p1 = [];
p2 = [];
p3 = [];
p4 = [];
p5 = [];
p6 = [];
p7 = [];
p8 = [];
p9 = [];
v1 = [];
v2 = [];
v3 = [];
v4 = [];
v5 = [];
v6 = [];
v7 = [];
v8 = [];
v9 = [];

% find opt solution...
A = [];
b = [];
Aeq = [];
beq = [];

% delta:
delta = [5*pi/180,... % q1 [rad]
        5*pi/180,... % q2
        10/1000,...  % R1 [m]
        10/1000,...  % R2
        10/1000,...  % R3
        10/1000,...  % Rc1
        10/1000,...  % Rc2
        10/1000,...  % Rc3
        10/1000,...  % e1
        10/1000,...  % e2
        50/1000,...  % e3
        50/1000].';  % e4

lb = x0 - delta; % lower bound
ub = x0 + delta; % upper bound

[x_sol,fval,exitflag,output,lambda,grad,hessian] = fmincon(@func_cost,x0,A,b,Aeq,beq,lb,ub,@func_nonlcon,options);

% save solution (& plot the mechanism) if converged
if exitflag == 1 || exitflag == 2
    
    % save
    dir = 'optsols\kin\';
    td = datetime;
    td_str = [num2str(td.Year),'_',num2str(td.Month),'_',num2str(td.Day),'_',num2str(td.Hour),'_',num2str(td.Minute),'_',num2str(td.Second)];
    save([dir,td_str,'.mat'],'x_sol','fval','exitflag','output','lambda','grad','hessian');
% %     
    % solve for unknowns
    % alpha0:
    %   alpha1
    %   alpha2
    %   alpha3
    %   alpha4
    %   alpha5
    %   alpha6
    % x: 
    %   q1
    %   q2
    %   R1 
    %   R2 
    %   R3 
    %   Rc1 
    %   Rc2 
    %   Rc3 
    %   e1 
    %   e2 
    %   e3 
    %   e4
    %
    % alpha_sol:
    %   alpha1
    %   alpha2
    %   alpha3
    %   alpha4
    %   alpha5
    %   alpha6
    [a_sol]=func_sol_mechanism(alpha0, x_sol(1), x_sol(2), x_sol(3), x_sol(4), x_sol(5), x_sol(6), x_sol(7), x_sol(8), x_sol(9), x_sol(10), x_sol(11), x_sol(12));
    
    % computing fwd kinematics terms
    % alpha_sol:
    %   alpha1
    %   alpha2
    %   alpha3
    %   alpha4
    %   alpha5
    %   alpha6
    % x: 
    %   q1
    %   q2
    %   R1 
    %   R2 
    %   R3 
    %   Rc1 
    %   Rc2 
    %   Rc3 
    %   e1 
    %   e2 
    %   e3 
    %   e4
    %
    % pi: joint positions
    % vi: joint velocities
    [p1,p2,p3,p4,p5,p6,p7,p8,p9,v1,v2,v3,v4,v5,v6,v7,v8,v9] = ...
        func_leg_fwd_kin(a_sol, x_sol(1), x_sol(2), x_sol(3), x_sol(4), x_sol(5), x_sol(6), x_sol(7), x_sol(8), x_sol(9), x_sol(10), x_sol(11), x_sol(12));
    
   
    % computing the jacobians 
    % alpha_sol:
    %   alpha1
    %   alpha2
    %   alpha3
    %   alpha4
    %   alpha5
    %   alpha6
    % x: 
    %   q1
    %   q2
    %   R1 
    %   R2 
    %   R3 
    %   Rc1 
    %   Rc2 
    %   Rc3 
    %   e1 
    %   e2 
    %   e3 
    %   e4
    %
    % J1: jacobian(p9,[q1;q2])
    % J2: jacobian(p9,[alpha1;alpha2;alpha3;alpha4;alpha5;alpha6])
    % Jc1: jacobian("holonomic constraints vector",[q1;q2])
    % Jc2: jacobian("holonomic constraints vector",[alpha1;alpha2;alpha3;alpha4;alpha5;alpha6])
    [J1_tmp,J2_tmp,Jc1_tmp,Jc2_tmp] = func_compute_jacs(a_sol, x_sol(1), x_sol(2), x_sol(3), x_sol(4), x_sol(5), x_sol(6), x_sol(7), x_sol(8), x_sol(9), x_sol(10), x_sol(11), x_sol(12));

    % compute torques
    A = J1_tmp - (J2_tmp/Jc2_tmp)*Jc1_tmp;

    % tau = [tau1, tau2]
    tau_tmp = -GRF.'*A;
    tau1 = tau_tmp(1);
    tau2 = tau_tmp(2);
    disp(['tau1 = ',num2str(tau1), ' N.m']);
    disp(['tau2 = ',num2str(tau2), ' N.m']);

    % plot leg mechanism
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
    
    % virtual leg
    line(ah,[(p1(1)+p2(1))/2,p9(1)],[(p1(2)+p2(2))/2,p9(2)],'Color','r','LineStyle','-.'); 
    
    
end





    % cost function
    %   x:
    %       q1
    %       q2
    %       R1
    %       R2
    %       R3
    %       Rc1
    %       Rc2
    %       Rc3
    %       e1
    %       e2
    %       e3
    %       e4
    %   f:
    %       tau.'*tau (minimize torque)
    function [f] = func_cost(x)
        
        % solve for unknowns
        % alpha0:
        %   alpha1
        %   alpha2
        %   alpha3
        %   alpha4
        %   alpha5
        %   alpha6
        % x: 
        %   q1
        %   q2
        %   R1 
        %   R2 
        %   R3 
        %   Rc1 
        %   Rc2 
        %   Rc3 
        %   e1 
        %   e2 
        %   e3 
        %   e4
        %
        % alpha_sol:
        %   alpha1
        %   alpha2
        %   alpha3
        %   alpha4
        %   alpha5
        %   alpha6
        [alpha_sol]=func_sol_mechanism(alpha0, x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10), x(11), x(12));
       
        % computing fwd kinematics terms
        % alpha_sol:
        %   alpha1
        %   alpha2
        %   alpha3
        %   alpha4
        %   alpha5
        %   alpha6
        % x: 
        %   q1
        %   q2
        %   R1 
        %   R2 
        %   R3 
        %   Rc1 
        %   Rc2 
        %   Rc3 
        %   e1 
        %   e2 
        %   e3 
        %   e4
        %
        % pi: joint positions
        % vi: joint velocities
        [p1,p2,p3,p4,p5,p6,p7,p8,p9,v1,v2,v3,v4,v5,v6,v7,v8,v9] = ...
            func_leg_fwd_kin(alpha_sol, x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10), x(11), x(12));
        
        % computing the jacobians 
        % alpha_sol:
        %   alpha1
        %   alpha2
        %   alpha3
        %   alpha4
        %   alpha5
        %   alpha6
        % x: 
        %   q1
        %   q2
        %   R1 
        %   R2 
        %   R3 
        %   Rc1 
        %   Rc2 
        %   Rc3 
        %   e1 
        %   e2 
        %   e3 
        %   e4
        %
        % J1: jacobian(p9,[q1;q2])
        % J2: jacobian(p9,[alpha1;alpha2;alpha3;alpha4;alpha5;alpha6])
        % Jc1: jacobian("holonomic constraints vector",[q1;q2])
        % Jc2: jacobian("holonomic constraints vector",[alpha1;alpha2;alpha3;alpha4;alpha5;alpha6])
        [J1,J2,Jc1,Jc2] = func_compute_jacs(alpha_sol, x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10), x(11), x(12));
        
        % compute torques
        tmp = J1 - (J2/Jc2)*Jc1;
        
        % tau = [tau1, tau2]
        tau = -GRF.'*tmp;
        
        % compute cost
        f = tau*tau.';
    end

    % nonlinear constraints
    %   x:
    %       q1
    %       q2
    %       R1
    %       R2
    %       R3
    %       Rc1
    %       Rc2
    %       Rc3
    %       e1
    %       e2
    %       e3
    %       e4
    %
    %   c: inequality constraints
    %   ceq: equality constraints
    function [c,ceq] = func_nonlcon(x)
        c = [];
        ceq = [];
        
        % compute nonlinear inequalities @ x
        c = [c,...
            p9(2) + 0.5,... % foot to hip distance
            x(1) - 160*pi/180,... % q1 < c
            100*pi/180 - x(1),... % c < q1
            x(2) - 160*pi/180,... % q2 < c
            100*pi/180 - x(2),... % c < q2
            60*pi/180 - alpha_sol(4),... % c < alpha4 
            0.25 - x(4),... % c < R2
            0.2 - x(7),... % c < Rc2
            0.07 - x(9),... % c < e1
            x(10) - 0.1,... % e2 < c
            0.07 - x(3),... % c < R1
            x(12) - 0.02 - x(11),... % e4 - c < e3
%             0.07 - x(6),... % c < Rc1
            ].';

        % compute nonlinear equalities @ x
        ceq = [ceq, (p1(1) + p2(1))/2 - p9(1),...
            x(6) - x(3),... % R1 = Rc1
            0.05 - x(12),... % e4 = c
%             0.05 - x(3),... % R1 = c
%             0.05 - x(6),... % Rc1 = c
            ].';
        
    end


end