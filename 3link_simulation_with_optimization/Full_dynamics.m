%Simulates full hybrid dynamics with parametrized output function
function Full_dynamics()
%Output of fmincon [q1, q1_dot, alpha]
f = [-0.2618 -2 0 0.20 0.5236 2.0944 2.50 2.618];
%fast dynamics
%f = [-0.4821   -6.3    0.1918    0.3672    0.6869  2.3175    2.2645    2.8316];
f = [-0.4431   -2.8   0.1702    0.5189    0.6723    1.9   2.5944    2.8558];

M = 4; delq = -deg2rad(30); x_plus = 0; 

global u2 time2
u2 = [];
time2 = [];

z_minus = f(1:2);
%alpha is the Bezier coeffs for q2, gamma for q3
alpha = [-f(5), -f(4), f(3:5)];
gamma = [-f(8)+2*f(6), -f(7)+2*f(6), f(6:8)];

a = [alpha, gamma];

%Mapping from z to x on boundary pg 140, 141
q2_minus = alpha(5);            %q2- = alpha(M)
q3_minus = gamma(5);            %q3- = gamma(M)

%d_dot- = M*(alpha(M) - alpha(M-1))*theta_dot-/(delta_theta)
dq2_minus = M*(alpha(5)-alpha(4))*z_minus(2)/delq;
dq3_minus = M*(gamma(5)-gamma(4))*z_minus(2)/delq;

%I.C. [thetas; velocities]
x_minus = [z_minus(1), q2_minus, q3_minus, z_minus(2), dq2_minus, dq3_minus];
x_plus = impact_map(x_minus);
x_plus = x_plus(1:6)';

delq = x_minus(1) - x_plus(1);

x_plus = [x_plus(1), bezier(0,4,alpha), bezier(0,4,gamma), x_plus(4), d_ds_bezier(0,4,alpha)*(x_plus(4)/delq), d_ds_bezier(0,4,gamma)*(x_plus(4)/delq)]';

tstart = 0; tfinal = 50;                    %max time per swing

options = odeset('Events',@events);    

time_out = tstart; states = x_plus.'; foot = [];
time_y = [];    y = []; 
phase_portrait = [x_minus; x_plus']; phase_x_minus = x_minus; phase_x_plus = x_plus';
%Simlulation
for i = 1:10                 %max number of steps allowed (if allowed by tfinal)
    
    [t,x] = ode45(@(t,x) dx_vector_field(t,x,a), [tstart tfinal], x_plus, options);
    nt = length(t);
    
    tstart = t(nt);
    
    %Save data
    time_out = [time_out; t(1:nt)]; time_y = [time_y; t(1:nt)];
    states = [states;x(1:nt,:)];    
    foot = [foot; (x(1,1)+x(end,1))/2*ones(length(x(:,1)),1) + states(end,1)  - x(end,1)]; %remember foot location
    %for error
    y = [y; output(x(1:nt,:),a,x_plus,delq)];
    %for phase portrait
    phase_x_minus = [phase_x_minus; x(nt,:)];
    temp_x = impact_map(phase_x_minus(end,:)); 
    phase_x_plus = [phase_x_plus; temp_x(1:6)];
    
    %Setting the new initial conditions based on impact map
    [x_plus,~] = impact_map(x(end,:));
    x_plus = x_plus(1:6);       % Only positions and velocities needed as initial conditions
        
end

[rows,~] = size(phase_x_minus);
for count = 2:rows
    phase_portrait = [phase_portrait; phase_x_minus(count,:); phase_x_plus(count,:)];
end
size(t)
%
%Plots phase portrait with impact
figure(1)
plot(states(:,1),states(:,4),'b-.')
hold on
[rows,~] = size(phase_portrait); 
for count = 1:2:rows
    plot(phase_portrait(count:count+1,1),phase_portrait(count:count+1,4),'rx-')
    if count == rows
        plot(phase_portrait(count,1),phase_portrait(count,4),'rx')
    end
end

hold off
legend('swing phase','impact','location','best','Interpreter','latex')
xlabel('q_1'); ylabel('$\dot{q_1}$','Interpreter','latex')
%}

%
%Plots of states
figure(2)
subplot(2,1,1)
plot(time_out,states(:,1),time_out,states(:,2),'-.',time_out,states(:,3),'--')
legend('\theta_1','\theta_2','\theta_3','location','best')
title('Joint angles')
subplot(2,1,2)
plot(time_out,states(:,4),time_out,states(:,5),'-.',time_out,states(:,6),'--')
legend('$\dot{\theta_1}$','$\dot{\theta_2}$','$\dot{\theta_3}$','location','best','Interpreter','latex')
title('Joint velocities')
%}

%
%Plot of error
figure(3)
plot(time_y,y(:,1))
hold on
plot(time_y,y(:,2))
plot(time_y,y(:,3))
plot(time_y,y(:,4))
hold off
legend('y1','y2','dy1','dy2','location','best')
%}

figure(5)
plot(time2,u2(1,:))
hold on
plot(time2,u2(2,:))
hold off
legend('u1','u2','location','best'), xlabel('time (sec)')

%{
%Stick figure plot
[hip_posx, leg1, leg2, torso] = motion(time_out,states);
hip = [hip_posx, zeros(size(hip_posx))];
figure(4)
plot(leg1(:,1),leg1(:,2),'o'), hold on
plot(leg2(:,1),leg2(:,2),'x')
plot(torso(:,1),torso(:,2),'x'),
[n,~] = size(hip);
for i = 1:n
    line([hip(i,1),leg1(i,1)], [hip(i,2),leg1(i,2)],'color','[0 0.4470 0.7410]','LineStyle','-')
    line([hip(i,1),leg2(i,1)], [hip(i,2),leg2(i,2)],'color','[0.8500    0.3250    0.0980]','LineStyle','-.')
    line([hip(i,1),torso(i,1)], [hip(i,2),torso(i,2)],'color','[0.9290    0.6940    0.1250]','LineStyle','--')
end
hold off
%}

function out = output(x,a,x_plus,delq)
    % Bezier coefficients 
    a2 = a(1:5); a3 = a(6:10);
    [m,~] = size(x);
    for k = 1:m
        q1 = x(k,1); q2 = x(k,2); q3 = x(k,3);
        s = (q1 - x_plus(1))/delq;         
        b2 = bezier(s,4,a2);
        b3 = bezier(s,4,a3);
        
        h = [q2 - b2; q3 - b3];         %y = h(x) = Hq - hd
        [D,C,G,~] = D_C_G_matrix(x(k,:));
        Fx = [x(k,4:6)'; D\(-C*x(k,4:6)'-G)];
        dh_dx = [ 0, 1, 0, 0, 0, 0;...
            0, 0, 1, 0, 0, 0];
        dh_dx(1,1) = -d_ds_bezier(s,4,a2)/delq;
        dh_dx(2,1) = -d_ds_bezier(s,4,a3)/delq;
        Lfh = dh_dx*Fx;
        
        out(k,:) = [h; Lfh]';
    end
end

%Event function
    function [limits,isterminal,direction] = events(~,x)
        
        q1 = x(1);
        s = (q1 - x_plus(1))/delq;   %normalized cyclic coordinate
        
        limits(1) = s-1;        %check when gait is at the end
        limits(2) = s+0.1;      %check if gait rolls back to 0
        isterminal = [1 1];    	% Halt integation
        direction = [];       	%The zero can be approached from either direction
    end

    function dx = dx_vector_field(time_in,x,a)
        % Computes vector field x_dot = f(x) + g(x)*u
        % Inputs:
        %       t: time
        %       x: all states [q, dq]
        %       a: bezier coefficient for q2 (1:5) and q3 (6:10)
        %
        % Outputs:
        %       dx = [dq, ddq]
        
        q1 = x(1); q2 = x(2); q3 = x(3); dq1 = x(4);
        
        [D,C,G,B] = D_C_G_matrix(x);
        
        Fx = [x(4:6); D\(-C*x(4:6)-G)];
        Gx = [zeros(3,2);D\B];
        
        % Bezier coefficients for q2
        a21 = a(1); a22 = a(2); a23 = a(3); a24 = a(4); a25 = a(5);
        a2 = a(1:5);
        % Bezier coefficients for q3
        a31 = a(6); a32 = a(7); a33 = a(8); a34 = a(9); a35 = a(10);
        a3 = a(6:10);
        
        %s variable used in bezier polynomial
        s = (q1 - x_plus(1))/delq;         %normalized cyclic coordinate
        
        b2 = bezier(s,4,a2);
        b3 = bezier(s,4,a3);
        
        %Output variable
        h = [q2 - b2; q3 - b3];         %y = h(x) = Hq - hd
        
        dh_dx = [ 0, 1, 0, 0, 0, 0;...
            0, 0, 1, 0, 0, 0];
        dh_dx(1,1) = -d_ds_bezier(s,4,a2)/delq;   %- db/ds*ds/dq1
        dh_dx(2,1) = -d_ds_bezier(s,4,a3)/delq;   
            
        Lfh = dh_dx*Fx;
        
        dLfh = [-dq1/delq^2*(12*a23*s^2 - 24*a24*s^2 + 12*a25*s^2 + 12*a21*(s - 1)^2 - 24*a22*(s - 1)^2 + 12*a23*(s - 1)^2 - 24*a24*s*(s - 1) - 12*a22*s*(2*s - 2) + 24*a23*s*(2*s - 2)), 0, 0,...
            -d_ds_bezier(s,4,a2)/delq, 1, 0;...
            -dq1/delq^2*(12*a33*s^2 - 24*a34*s^2 + 12*a35*s^2 + 12*a31*(s - 1)^2 - 24*a32*(s - 1)^2 + 12*a33*(s - 1)^2 - 24*a34*s*(s - 1) - 12*a32*s*(2*s - 2) + 24*a33*s*(2*s - 2)), 0, 0,...
            -d_ds_bezier(s,4,a3)/delq, 0, 1];      
        
        %
        epsilon = 0.05; alp = 0.5;
        
        %scaling
        Lfh = epsilon*Lfh;
        
        phi(1) = h(1) + 1/(2 - alp)*sign(Lfh(1))*abs(Lfh(1))^(2-alp);
        phi(2) = h(2) + 1/(2 - alp)*sign(Lfh(2))*abs(Lfh(2))^(2-alp);
        
        psi(1,1) = -sign(Lfh(1))*abs(Lfh(1))^alp - sign(phi(1))*abs(phi(1))^(alp/(2-alp));
        psi(2,1) = -sign(Lfh(2))*abs(Lfh(2))^alp - sign(phi(2))*abs(phi(2))^(alp/(2-alp));
        
        v = 1/epsilon^2*psi;
        %}
        %{
        %PD control
        eps = 0.05;
        Kp = diag([10/eps,10/eps]);
        Kd = diag([10/eps^2,10/eps^2]);
        v = -(Kp*h + Kd*Lfh);
        %}
        
        u = (dLfh*Gx)\(v - dLfh*Fx);    
        dx = Fx + Gx*u;
        u2 = [u2, u];
        time2 = [time2, time_in];
    end
    

    function [hip_posx, leg1, leg2, torso] = motion(t,x)
        % Estimate hip horizontal position by estimating integral of hip velocity
        
        [r,~,~,~,l,~] = model_params_3link;
        
        hip_velx = r*cos(x(:,1)).*x(:,4);     %r*cos(q1)*dq1
        
        [n,~]=size(x);
        hip_posx = zeros(n,1);
        
        for j=2:n
            hip_posx(j)=hip_posx(j-1)+(t(j)-t(j-1))*hip_velx(j-1,1);
        end
        
        leg1 = [hip_posx + r*sin(x(:,1)), -r*cos(x(:,1))];
        leg2 = [hip_posx - r*sin(x(:,2) + x(:,1)), -r*cos(x(:,2) + x(:,1))];
        torso = [hip_posx + l*sin(pi - x(:,3) - x(:,1)), l*cos(pi - x(:,3) - x(:,1))];
        
    end
end