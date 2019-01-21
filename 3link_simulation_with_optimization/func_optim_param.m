% Optimize intial conditions as ZD states and Bezier coeffs using fmincon
function [f,J] = func_optim_param()
z_minus = [-0.2618 -2];         %[q1, dq1]
%J = 0;                  %initial cost
z_sol = [];
delq = deg2rad(30);
%Bezier coefficients
alpha = [0 0.20 0.5236];            %for q2 - alpha 3 - 5
gamma = [2.0944 2.50 2.618];        %for q3 - alpha 3 - 5

y0 = [z_minus, alpha, gamma];   %parameters that need to be optimized

%boundary constraints starting from [q1 q1_dot, alpha q2, alpha q3]
eps = 0.25;
lb = [max(-deg2rad(30),z_minus(1)-eps), max(-10,z_minus(2)-3),...
    alpha(1)-eps, max(0,alpha(2)-eps), max(0,alpha(3)-eps),...
    max(0,gamma(1)-eps), max(0,gamma(2)-eps), max(0,gamma(3)-eps)];
ub = [0, -1,...
    alpha(1)+eps, min(deg2rad(45),alpha(2)+eps), min(deg2rad(45),alpha(3)+eps),...
    min(pi,gamma(1)+eps), min(pi,gamma(2)+eps), min(pi,gamma(3)+eps)];

%using only bounds and nonlinear contraints
opts = optimoptions('fmincon','Algorithm','interior-point');
[f, J] = fmincon(@(y)cost_ZD(y),y0,[],[],[],[],lb,ub,@confuneq,opts);

    function [c,ceq] = confuneq(y)
        
        z_minus = y(1:2);       %Pre-impact states
        z_f = z_sol(end,:);     %final values after swing phase and right before the next impact
        
        % Nonlinear inequality constraints
        c = [];
        % Nonlinear equality constraints
        ceq = [norm(z_minus(1)-z_f(1)); norm(z_minus(2)-z_f(2))];
        
    end

    function J = cost_ZD(y)
        
        J = 0;
        
        z_minus = y(1:2);
        %alpha is the Bezier coeffs for q2, gamma for q3
        alpha = y(3:5);         %alpha 3 - 5 for q2
        gamma = y(6:8);         %alpha 3 - 5 for q3
        
        M = 4;
        delq = 2*abs(z_minus(1));              %difference between min and max q1
        %Mapping from z to x on boundary pg 140, 141
        q2_minus = alpha(3);            %q- = alpha(M) (end of gait)
        q3_minus = gamma(3);
        %d_dot- = M*(alpha(M) - alpha(M-1))*theta_dot-/(delta_theta)
        dq2_minus = M*(alpha(3)-alpha(2))*z_minus(2)/delq;
        dq3_minus = M*(gamma(3)-gamma(2))*z_minus(2)/delq;
        
        x_minus = [z_minus(1), q2_minus, q3_minus, z_minus(2), dq2_minus, dq3_minus];
        
        %Applying impact map
        [x_plus, ~] = impact_map(x_minus);
        %Inverse map to get ZD states
        z_plus = [x_plus(1), x_plus(4)];
        dq2 = x_plus(5); dq3 =  x_plus(6);
        
        % Need the other bezier coefficients now
        
        a0 = -alpha(3); g0 = -gamma(3)+2*gamma(1);   %a(0) = -a(5); g(0) = -g(4)+2*g(2)
        %Enforcing slope for the 1st 2 coefficients to be equal to qb_dot:
        %M*(alpha(1) - alpha(0))*theta_dot-/(delta_theta) = qb_dot, so:
        a1 = -alpha(2);   % -a5;
        g1 = -gamma(3)+2*gamma(1);   % g(1)-g(0) = g(4)-g(3);
        
        a = [a0, a1, alpha, g0, g1, gamma];
        
        % compute delq
        delq = z_plus(1)-z_minus(1);
        
        tstart = 0; tfinal = 1;                    %max time per swing
        
        refine = 4; options = odeset('Events',@events,'Refine',refine);
        
        [~,z_sol] = ode45(@(t,z) ZD_states(t,z,a), [tstart tfinal], z_plus, options);
        
        % J = z_sol(end,3);       %last entry of the control input norm integral? Should it not be the last entry?
        
        %{
%Alternate way to get cost, by summing control input vector, need a
%function that computes control input
[m1,~] = size(z);       %get total entries for a single sim
%map ZD coordinates vector z to vector full_x
for j = 1:m1
    full_x(j,:) = map_z_to_x(z,a);
end
[m2,~] = size(full_x);  %m2 is same as m1 without sampling
J = 0;
%get u for each row of x, get norm 2 of row, then sum to get J
for i = 1:m2
    u(i,:) = control_input(full_x(i,:), a);
    normed_u(i,:) = norm(u(i,:));
    J = J + normed_u(i,:);
end
        %}
        
        %
        %Event function
        function [limits,isterminal,direction] = events(~,z)
            
            q1 = z(1);
            s = (z_plus(1) - q1)/delq;   %normalized general coordinate
            
            limits(1) = s-1; 	%check when gait is at the end
            limits(2) = s;      %check if gait rolls back to 0
            isterminal = [1 1];    	% Halt integation
            direction = [];       	%The zero can be approached from either direction
            
        end
        
        %Zero Dynamics
        function dz = ZD_states(~,z,a)
            
            x = map_z_to_x(z,a);
            [D,C,G,B] = state_matrix(x);        %Get state matrix using full states x
            
            % Bezier coefficients for q2
            a21 = a(1); a22 = a(2); a23 = a(3); a24 = a(4); a25 = a(5);
            a_2 = a(1:5);
            % Bezier coefficients for q3
            a31 = a(6); a32 = a(7); a33 = a(8); a34 = a(9); a35 = a(10);
            a_3 = a(6:10);
            
            q1 = z(1);
            dq1 = z(2);
            
            D1 = D(1, 1);
            D2 = D(1, 2:3);
            D3 = D(2:3, 1);
            D4 = D(2:3, 2:3);
            
            H1 = C(1,1)*dq1 + G(1,1);
            H2 = C(2:3,2:3)*x(5:6)' + [G(2,1); G(3,1)];
            
            %dq1,q1
            s = (z_plus(1)-q1)/delq;   %normalized general coordinate
            %x(2),s
            % d/ds(db/ds*s_dot)
            dLsb2 = -dq1/delq*(3*s^2*(4*a24 - 4*a25) - s^2*(12*a23 - 12*a24) - 3*(s - 1)^2*(4*a21 - 4*a22) +...
                (s - 1)^2*(12*a22 - 12*a23) - 2*s*(s - 1)*(12*a23 - 12*a24) + s*(2*s - 2)*(12*a22 - 12*a23));
            
            dLsb3 = -dq1/delq*(3*s^2*(4*a34 - 4*a35) - s^2*(12*a33 - 12*a34) - 3*(s - 1)^2*(4*a31 - 4*a32) +...
                (s - 1)^2*(12*a32 - 12*a33) - 2*s*(s - 1)*(12*a33 - 12*a34) + s*(2*s - 2)*(12*a32 - 12*a33));
            
            % beta1 =  d/ds(db/ds*s_dot)*s_dot
            beta1 = [dLsb2; dLsb3]*dq1/delq;
            
            % beta2 =  db/ds*s_dot
            db_ds2 = d_ds_bezier(s,4,a_2);
            db_ds3 = d_ds_bezier(s,4,a_3);
            
            beta2 = [db_ds2; db_ds3]/delq;
            
            ddq1 = (D1 + D2*beta2)\(-D2*beta1 - H1);
            
            u = B(2:3,1:2)\((D3 + D4*beta2)*ddq1 + (D4*beta1 + H2));
            
            dz(1) = z(2);
            dz(2) = (D1 + D2*beta2)\(-D2*beta1 - H1);
            
            %Setting the derivative of the cost function as a state, from pg155
            %dz(3) = norm(u);
            J = J + norm(u);
            dz = [dz(1), dz(2)]';
            
        end
        %}
        
        function q = map_z_to_x(z,a)
            
            q1 = z(1);
            dq1 = z(2);
            a2 = a(1:5); a3 = a(6:end);
            
            M = 4;
            
            s = (z_plus(1) - q1)/delq;   %normalized general coordinate
            
            q2 = bezier(s,M,a2);
            q3 = bezier(s,M,a3)';
            
            dq2 = d_ds_bezier(s,M,a2)*dq1/delq;
            dq3 = d_ds_bezier(s,M,a3)*dq1/delq;
            
            q = [z(1), q2, q3, z(2), dq2, dq3];
            
        end
        
    end

end