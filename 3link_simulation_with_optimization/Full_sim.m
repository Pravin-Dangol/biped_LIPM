%Simulates full hybrid dynamics with parametrized output function
function Full_sim()
%Output of fmincon [q1, q1_dot, alpha]
f = [-0.2618 -4 0 0.20 0.5236 2.0944 2.50 2.618];
f = [-0.4821   -2.197    0.1918    0.3672    0.6869  2.3175    2.2645    2.8316];
M = 4; delq = deg2rad(30); x_plus = 0; y = []; dy = [];

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
%x_minus = map_z_to_x([x_minus(1),x_minus(4)],a);
%Apply impact map
x_plus = impact_map(x_minus);
x_plus = x_plus(1:6)';


delq = x_plus(1) - x_minus(1);

tstart = 0; tfinal = 50;                    %max time per swing

refine = 4; options = odeset('Events',@events,'Refine',refine);    %'OutputFcn',@odeplot,'OutputSel',1,

time_out = tstart; states = x_plus.'; foot = [];
time_y = [];
%Simlulation
for i = 1:2                 %max number of steps allowed (if allowed by tfinal)
    
    [t,x] = ode45(@(t,x) dx_vector_field(t,x,a), [tstart tfinal], x_plus, options);
    nt = length(t);
    
    %Setting the new initial conditions based on impact map
    [x_plus,~] = impact_map(x(end,:));
    x_plus = x_plus(1:6);       % Only positions and velocities needed as initial conditions
    
    %
    %Stop conditions
    if tstart >= tfinal || nt-refine < 1
        break
    end
    %}
    
    options = odeset(options,'InitialStep',t(nt)-t(nt-refine),'MaxStep',t(nt)-t(1));
    tstart = t(nt);
    
    %Save data
    time_out = [time_out; t(2:nt)]; time_y = [time_y; t(1:nt)];
    states = [states;x(2:nt,:)];    
    foot = [foot; (x(1,1)+x(end,1))/2*ones(length(x(:,1)),1) + states(end,1)  - x(end,1)]; %remember foot location
    y = [y; output(x(1:nt,:),a,x_plus,delq)];
    
end

%{
%Plots
figure(1)
subplot(2,1,1)
plot(time_out,states(:,1),time_out,states(:,2),'-.',time_out,states(:,3),'--')
legend('\theta_1','\theta_2','\theta_3')
title('Joint angles')
subplot(2,1,2)
plot(time_out,states(:,4),time_out,states(:,5),'-.',time_out,states(:,6),'--')
legend('$\dot{\theta_1}$','$\dot{\theta_2}$','$\dot{\theta_3}$','Interpreter','latex')
title('Joint velocities')
%}



%
plot(time_y,y(:,1))
hold on
plot(time_y,y(:,2))
plot(time_y,y(:,3))
plot(time_y,y(:,4))
hold off
legend('y1','y2','dy1','dy2')
%}

%{
%animation
[hip_posx, leg1, leg2, torso] = motion(time_out,states);
hip = [hip_posx, zeros(size(hip_posx))];
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
    a2 = a(1:5); a3 = a(6:10);
    [m,~] = size(x);
    for k = 1:m
        q1 = x(k,1); q2 = x(k,2); q3 = x(k,3);
        s = (x_plus(1) - q1)/delq;         %normalized cyclic coordinate
        b2 = bezier(s,4,a2);
        b3 = bezier(s,4,a3);
        
        h = [q2 - b2; q3 - b3];         %y = h(x) = Hq - hd
        [D,C,G,~] = state_matrix(x(k,:));
        Fx = [x(k,4:6)'; D\(-C*x(k,4:6)'-G)];
        dh_dx = [ 0, 1, 0, 0, 0, 0;...
            0, 0, 1, 0, 0, 0];
        dh_dx(1,1) = -d_ds_bezier(s,4,a2);
        dh_dx(2,1) = -d_ds_bezier(s,4,a3);
        Lfh = dh_dx*Fx;
        
        out(k,:) = [h; Lfh]';
    end
end

%Event function
    function [limits,isterminal,direction] = events(~,x)
        
        q1 = x(1);
        s = (x_plus(1) - q1)/delq;   %normalized cyclic coordinate
        
        limits(1) = s-1;        %check when gait is at the end
        limits(2) = s+0.1;      %check if gait rolls back to 0
        isterminal = [1 1];    	% Halt integation
        direction = [];       	%The zero can be approached from either direction
    end

    function dx = dx_vector_field(~,x,a)
        
        %Computes vector field x_dot = f(x) + g(x)*u
        %Required inputs: x - all states [q; q_dot] and a - Bezier coeffs
        
        t_temp = num2cell(x(1:3)); [q1, q2, q3] = t_temp{:};
        dt_temp = num2cell(x(4:6)); [dq1, ~, ~] = dt_temp{:};
        
        [D,C,G,B] = state_matrix(x);
        
        Fx = [x(4:6); D\(-C*x(4:6)-G)];
        Gx = [zeros(3,2);D\B];
        
        % Bezier coefficients for q2
        a21 = a(1); a22 = a(2); a23 = a(3); a24 = a(4); a25 = a(5);
        a2 = a(1:5);
        % Bezier coefficients for q3
        a31 = a(6); a32 = a(7); a33 = a(8); a34 = a(9); a35 = a(10);
        a3 = a(6:10);
        
        %s variable used in bezier polynomial
        s = (x_plus(1) - q1)/delq;         %normalized cyclic coordinate
        
        b2 = bezier(s,4,a2);
        b3 = bezier(s,4,a3);
        
        %Output variable
        h = [q2 - b2; q3 - b3];         %y = h(x) = Hq - hd
        
        dh_dx = [ 0, 1, 0, 0, 0, 0;...
            0, 0, 1, 0, 0, 0];
        dh_dx(1,1) = -d_ds_bezier(s,4,a2);
        dh_dx(2,1) = -d_ds_bezier(s,4,a3);
        
        Lfh = dh_dx*Fx;
        
        dLfh = [-dq1*(12*a23*s^2 - 24*a24*s^2 + 12*a25*s^2 + 12*a21*(s - 1)^2 - 24*a22*(s - 1)^2 + 12*a23*(s - 1)^2 - 24*a24*s*(s - 1) - 12*a22*s*(2*s - 2) + 24*a23*s*(2*s - 2)), 0, 0, -d_ds_bezier(s,4,a2), 1, 0;...
            -dq1*(12*a33*s^2 - 24*a34*s^2 + 12*a35*s^2 + 12*a31*(s - 1)^2 - 24*a32*(s - 1)^2 + 12*a33*(s - 1)^2 - 24*a34*s*(s - 1) - 12*a32*s*(2*s - 2) + 24*a33*s*(2*s - 2)), 0, 0, -d_ds_bezier(s,4,a3), 0, 1];
        %{
        epsilon = 0.1; alp = 0.9;
        
        %scaling
        Lfh = epsilon*Lfh;
        
        phi(1) = h(1) + 1/(2 - alp)*sign(Lfh(1))*abs(Lfh(1))^(2-alp);
        phi(2) = h(2) + 1/(2 - alp)*sign(Lfh(2))*abs(Lfh(2))^(2-alp);
        
        psi(1,1) = -sign(Lfh(1))*abs(Lfh(1))^alp - sign(phi(1))*abs(phi(1))^(alp/(2-alp));
        psi(2,1) = -sign(Lfh(2))*abs(Lfh(2))^alp - sign(phi(2))*abs(phi(2))^(alp/(2-alp));
        
        v = 1/epsilon^2*psi;
        %}
        %
        %PD control
        eps=0.01;
        Kp = diag([10/eps,1/eps]);
        Kd = diag([1/eps^2,1/eps^2]);
        v = -(Kp*h + Kd*dh_dx*Fx);
        %}
        
        %{
        % Bezier coefficients for q2
        D1 = D(1, 1);
        D2 = D(1, 2:3);
        D3 = D(2:3, 1);
        D4 = D(2:3, 2:3);
        
        H1 = C(1,1)*dq1 + G(1,1);
        H2 = C(2:3,2:3)*x(5:6) + [G(2,1); G(3,1)];
        
        dLsb2 = -dq1/delq*(3*s^2*(4*a24 - 4*a25) - s^2*(12*a23 - 12*a24) - 3*(s - 1)^2*(4*a21 - 4*a22) +...
            (s - 1)^2*(12*a22 - 12*a23) - 2*s*(s - 1)*(12*a23 - 12*a24) + s*(2*s - 2)*(12*a22 - 12*a23));
        
        dLsb3 = -dq1/delq*(3*s^2*(4*a34 - 4*a35) - s^2*(12*a33 - 12*a34) - 3*(s - 1)^2*(4*a31 - 4*a32) +...
            (s - 1)^2*(12*a32 - 12*a33) - 2*s*(s - 1)*(12*a33 - 12*a34) + s*(2*s - 2)*(12*a32 - 12*a33));
        
        % beta1 =  d/ds(db/ds*s_dot)*s_dot
        beta1 = [dLsb2; dLsb3]*dq1/delq;
        
        % beta2 =  db/ds*s_dot
        db_ds2 = d_ds_bezier(s,4,a2);
        db_ds3 = d_ds_bezier(s,4,a3);
        
        beta2 = [db_ds2; db_ds3]/delq;
        
        ddq1 = (D1 + D2*beta2)\(-D2*beta1 - H1);
        
        u = B(2:3,1:2)\((D3 + D4*beta2)*ddq1 + (D4*beta1 + H2));
        %}
        
        u = (dLfh*Gx)\(v - dLfh*Fx);    %u = LgLfh^-1*(v - L2fh) 
        
        dx = Fx + Gx*u;
        
    end
    
%
    function q = map_z_to_x(z,a)
        
        q1 = z(1);
        dq1 = z(2);
        a2 = a(1:5); a3 = a(6:end);
        
        M = 4;
        
        s = (x_plus(1) - q1)/delq;   %normalized general coordinate
        
        q2 = bezier(s,M,a2);
        q3 = bezier(s,M,a3)';
        
        dq2 = d_ds_bezier(s,M,a2)*dq1/delq;
        dq3 = d_ds_bezier(s,M,a3)*dq1/delq;
        
        q = [z(1), q2, q3, z(2), dq2, dq3];
        
    end
%}

    function [hip_posx, leg1, leg2, torso] = motion(t,x)
        
        [r,~,~,~,l,~] = model_params_3link;
        
        hip_velx = cos(x(:,1)).*x(:,4);
        
        [n,~]=size(x);
        hip_posx = zeros(n,1);
        % Estimate hip horizontal position by estimating integral of hip velocity
        for j=2:n
            hip_posx(j)=hip_posx(j-1)+(t(j)-t(j-1))*hip_velx(j-1,1);
        end
        
        leg1 = [hip_posx + r*sin(x(:,1)), -r*cos(x(:,1))];
        leg2 = [hip_posx - r*sin(x(:,2)), -r*cos(x(:,2))];
        torso = [hip_posx + l*sin(x(:,3)), l*cos(x(:,3))];
        
    end
end