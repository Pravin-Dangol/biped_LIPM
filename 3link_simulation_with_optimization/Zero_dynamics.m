%Simulation of a single step of ZD to compare with full dynamics
function Zero_dynamics()

f = [-0.2618 -2 0 0.20 0.5236 2.0944 2.50 2.618];
f = [-0.4821   -4.687    0.1918    0.3672    0.6869  2.3175    2.2645    2.8316];
%f = [-0.4431   -1.1779   0.1702    0.5189    0.6723    1.9   2.5944    2.8558];
delq = -2*abs(deg2rad(15)); 
z_plus = 0; M = 4; 

alpha = [-f(5), -f(4), f(3:5)];      
gamma = [-f(8)+2*f(6), -f(7)+2*f(6), f(6:8)];

a = [alpha, gamma];

z_minus = [f(1),f(2)];

%Mapping from z to x on boundary pg 140, 141
q2_minus = alpha(5);            %q2- = alpha(M)
q3_minus = gamma(5);            %q3- = gamma(M)
%d_dot- = M*(alpha(M) - alpha(M-1))*theta_dot-/(delta_theta)
dq2_minus = M*(alpha(5)-alpha(4))*z_minus(2)/delq;
dq3_minus = M*(gamma(5)-gamma(4))*z_minus(2)/delq;

x_minus = [z_minus(1), q2_minus, q3_minus, z_minus(2), dq2_minus, dq3_minus];

x_plus = impact_map(x_minus);
x_plus = x_plus(1:6);

z_plus = [x_plus(1), x_plus(4)];

delq = z_minus(1) - z_plus(1);

tstart = 0; tfinal = 50;                    %max time per swing
time_out = tstart; states = z_plus; 

impact_ = [z_minus; z_plus]; 
states_full = map_z_to_x(z_plus,a);

refine = 4; options = odeset('Events',@events,'Refine',refine);
for i = 1:1  
    
    [t,z] = ode45(@(t,z) ZD_states(t,z,a), [tstart tfinal], z_plus, options);
        
    nt = length(t); time_out = [time_out; t(2:nt)];
    states = [states;z(1:nt,:)];
    
    options = odeset(options,'InitialStep',t(nt)-t(nt-refine),'MaxStep',t(nt)-t(1));
    tstart = t(nt);
    
    x_minus = [x_minus; map_z_to_x(z(nt,:),a)];
    z_minus = [x_minus(end,1), x_minus(end,4)];
    temp_x = impact_map(x_minus(end,:)); 
    x_plus = [x_plus; temp_x(1:6)];
    z_plus = [x_plus(end,1), x_plus(end,4)];
    delq = z_minus(1) - z_plus(1);
    
    temp = zeros(1,6);
    for j = 2:nt
        temp(j-1,:) = map_z_to_x(z(j,:),a);
    end
    
    states_full = [states_full;temp];
         
end

[m,~] = size(x_minus);
for j = 2:m
    impact_ = [impact_; x_minus(j,1), x_minus(j,4); x_plus(j,1), x_plus(j,4)];
end

%
%Plots phase portrait with impact
figure(1)
plot(states(:,1),states(:,2),'b-.')
hold on
[m,~] = size(impact_); count = 1;
for k = 1:m-2
    plot(impact_(count:count+1,1),impact_(count:count+1,2),'rx-')
    count =  k+2;
end
hold off
legend('swing phase','impact','location','best','Interpreter','latex')
xlabel('q_1'); ylabel('$\dot{q_1}$','Interpreter','latex')
%}

%
%Plot all states
figure(2)
subplot(2,1,1)
plot(time_out,states_full(:,1),time_out,states_full(:,2),'-.',time_out,states_full(:,3),'--')
legend('q_1','q_2','q_3')
title('Joint angles')
subplot(2,1,2)
plot(time_out,states_full(:,4),time_out,states_full(:,5),'-.',time_out,states_full(:,6),'--')
legend('$\dot{q_1}$','$\dot{q_2}$','$\dot{q_3}$','location','best','Interpreter','latex')
title('Joint velocities')
%}

%Event function
    function [limits,isterminal,direction] = events(~,z)
        
        q1 = z(1);
        s = (q1 - z_plus(1))/delq;   %normalized general coordinate
        
        limits(1) = s-1; 	%check when gait is at the end
        limits(2) = s;      %check if gait rolls back to 0
        isterminal = [1 1];    	% Halt integation
        direction = [];       	%The zero can be approached from either direction
        
    end

%Zero Dynamics
    function dz = ZD_states(~,z,a)
        
        x = map_z_to_x(z,a);
        [D,C,G,~] = state_matrix(x);
        
        a21 = a(1); a22 = a(2); a23 = a(3); a24 = a(4); a25 = a(5); a_2 = a(1:5);
        a31 = a(6); a32 = a(7); a33 = a(8); a34 = a(9); a35 = a(10); a_3 = a(6:10);
        
        q1 = z(1); dq1 = z(2);
        
        D1 = D(1, 1);
        D2 = D(1, 2:3); %D3 = D(2:3, 1); %D4 = D(2:3, 2:3);
        
        H1 = C(1,1)*dq1 + G(1,1); %H2 = C(2:3,2:3)*x(5:6)' + [G(2,1); G(3,1)];
        
        s = (q1 - z_plus(1))/delq;   %normalized general coordinate

        dLsb2 = -dq1/delq*(3*s^2*(4*a24 - 4*a25) - s^2*(12*a23 - 12*a24) - 3*(s - 1)^2*(4*a21 - 4*a22) +...
            (s - 1)^2*(12*a22 - 12*a23) - 2*s*(s - 1)*(12*a23 - 12*a24) + s*(2*s - 2)*(12*a22 - 12*a23));
        
        dLsb3 = -dq1/delq*(3*s^2*(4*a34 - 4*a35) - s^2*(12*a33 - 12*a34) - 3*(s - 1)^2*(4*a31 - 4*a32) +...
            (s - 1)^2*(12*a32 - 12*a33) - 2*s*(s - 1)*(12*a33 - 12*a34) + s*(2*s - 2)*(12*a32 - 12*a33));
        
        beta1 = [dLsb2; dLsb3]*dq1/delq;
        
        db_ds2 = d_ds_bezier(s,4,a_2); db_ds3 = d_ds_bezier(s,4,a_3);
        
        beta2 = [db_ds2; db_ds3]/delq;
        
        dz(1) = z(2);
        dz(2) = (D1 + D2*beta2)\(-D2*beta1 - H1);
        %dz(2) = -G(1,1);
        
        dz = [dz(1), dz(2)]';
        
    end

    function q = map_z_to_x(z,a)
        
        q1 = z(1);
        dq1 = z(2);
        a2 = a(1:5); a3 = a(6:end);
        
        M = 4;
        
        s = (q1 - z_plus(1))/delq;   %normalized general coordinate
        
        q2 = bezier(s,M,a2);
        q3 = bezier(s,M,a3)';
        
        dq2 = d_ds_bezier(s,M,a2)*dq1/delq;
        dq3 = d_ds_bezier(s,M,a3)*dq1/delq;
        
        q = [z(1), q2, q3, z(2), dq2, dq3];
        
    end

end