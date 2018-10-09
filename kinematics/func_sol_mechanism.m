% solve for the unknowns in the leg mechanism
% 12/22/2017 Writen by Alirza Ramezani, Pasadena, CA
% x0:
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

% x_sol:
% alpha1
% alpha2
% alpha3
% alpha4
% alpha5
% alpha6

function [x_sol]=func_sol_mechanism(x0, q1, q2, R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4)


% solver settings...
options = optimoptions('fsolve',...
                            'diagnostics','on',...
                                   'display','iter',...
                                        'MaxFunctionEvaluations',5000,...
                                            'MaxIterations',400,...
                                                'OptimalityTolerance',1e-6,...
                                                    'UseParallel',false,...
                                                        'ScaleProblem','none'); 
           
% options = optimoptions('fsolve','Display','off',...
%                             'MaxFunctionEvaluations',400,...
%                                     'MaxIterations',400);

% solve for unknows in the first loop ...
% alpha1
% alpha2
% alpha3
disp('solving for unknows in the first loop ...');
[tmp1,~,exitflag] = fsolve(@func_constraints_first_loop,x0(1:3),options); % Call solver

% check the exitflag...
% 1:  equation solved. first-order optimality is small.
% 2:  equation solved. change in x smaller than the specified tolerance.
% 3:  equation solved. change in residual smaller than the specified tolerance.
% 4:  equation solved. magnitude of search direction smaller than specified tolerance.
% others, equation not solved.
if (exitflag==1)||(exitflag==2)||(exitflag==3)||(exitflag==4)
    disp('solution found.');
else
    disp('solution not found.');
%     break;
end

% alpha2 is required in the second loop...
ALPHA2 = tmp1(2);

% solve for unknows in the second loop ...
% alpha4
% alpha5
% alpha6
disp('solving for unknows in the second loop ...');
[tmp2,~,exitflag] = fsolve(@func_constraints_second_loop,x0(4:end),options); % Call solver

% check the exitflag...
% 1:  equation solved. first-order optimality is small.
% 2:  equation solved. change in x smaller than the specified tolerance.
% 3:  equation solved. change in residual smaller than the specified tolerance.
% 4:  equation solved. magnitude of search direction smaller than specified tolerance.
% others, equation not solved.
if (exitflag==1)||(exitflag==2)||(exitflag==3)||(exitflag==4)
    disp('solution found.');
else
    disp('solution not found.');
%     break;
end

x_sol = [tmp1.';tmp2.'];

% enforces the constraints in the first loop of LEO's leg mechanism...
%x:
% alpha1
% alpha2
% alpha3
    function [c]=func_constraints_first_loop(x)
        
        alpha1 = x(1);
        alpha2 = x(2);
        alpha3 = x(3);
        
        % compute constraints (symbolically generated ... look at "generate_constraints.m")
        c = zeros(3,1);
        c(1) = e2*cos(alpha1 + alpha2 - q2) - e1 + R1*cos(alpha1 + alpha2 + alpha3 - q2) + Rc1*cos(q2) + Rc2*cos(alpha1 - q2) - e3*cos(alpha1 - q2);
        c(2) = e2*sin(alpha1 + alpha2 - q2) + R1*sin(alpha1 + alpha2 + alpha3 - q2) - Rc1*sin(q2) + Rc2*sin(alpha1 - q2) - e3*sin(alpha1 - q2);
        c(3) = alpha1 + alpha2 + alpha3 - q1 - q2;
 
    end


% enforces the constraints in the second loop of LEO's leg mechanism...
%x:
% alpha4
% alpha5
% alpha6
    function [c]=func_constraints_second_loop(x)
        
        alpha4 = x(1);
        alpha5 = x(2);
        alpha6 = x(3);
        
        % alpha2 from solving the first loop ...
        alpha2 = ALPHA2;
        
        % compute constraints (symbolically generated ... look at "generate_constraints.m")
        c = zeros(3,1);
        c(1) = e2 - R2 - Rc3*cos(alpha4 + alpha5) - e4*cos(alpha4) - e3*cos(alpha4 + alpha5 + alpha6);
        c(2) = - Rc3*sin(alpha4 + alpha5) - e4*sin(alpha4) - e3*sin(alpha4 + alpha5 + alpha6);
        c(3) = alpha2 + alpha4 + alpha5 + alpha6 - 2*pi;
    end

end