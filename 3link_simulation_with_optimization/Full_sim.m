%Simulates full hybrid dynamics with parametrized output function
function Full_sim()
%Output of fmincon [q1, q1_dot, alpha,J]
f = [-0.2618 -1.2 0 0.20 0.5236 2.0944 2.50 2.618];
M = 4; delq = deg2rad(30);              %difference between min and max q1

z_minus = f(1:2);
%alpha is the Bezier coeffs for q2, gamma for q3
alpha = f(3:5);         %alpha 3 - 5 for q2
gamma = f(6:8);         %alpha 3 - 5 for q3

%Mapping from z to x on boundary pg 140, 141
%At the boundry (right before impact, q1_minus)
q2_minus = alpha(3);            %q- = alpha(M) (end of gait)
q3_minus = gamma(3);

a0 = -alpha(3); g0 = gamma(3)-1.0472;       

%qb_dot = db/ds*s_dot = M*(alpha(4) - alpha(3))*theta_dot-/(delta_theta)
dq2_minus = M*(alpha(3) - alpha(2))*z_minus(2)/delq;
dq3_minus = M*(gamma(3) - gamma(2))*z_minus(2)/delq;
%Enforcing slope for the 1st 2 coefficients to be equal to qb_dot:
%M*(alpha(1) - alpha(0))*theta_dot-/(delta_theta) = qb_dot, so:
a1 = dq2_minus*delq/(M*z_minus(2)) + a0;
g1 = dq3_minus*delq/(M*z_minus(2)) + g0;

alpha = [-f(5), -f(4), f(3:5)];      
gamma = [f(8)-f(6)/2, 2*f(8)-f(6)/2-f(7), f(6:8)];

a = [alpha, gamma];

%I.C. [thetas; velocities]
x0 = [z_minus(1), q2_minus, q3_minus, z_minus(2), dq2_minus, dq3_minus];

%Apply impact map
x_plus = impact_map(x0);
x_plus = x_plus(1:6)'; 

z_plus = [x_plus(1), x_plus(4)];
delq = z_plus(1) - z_minus(1);


tstart = 0; tfinal = 10;                    %max time per swing

refine = 4; options = odeset('Events',@events,'Refine',refine);    %'OutputFcn',@odeplot,'OutputSel',1,

time_out = tstart; states = x_plus.'; foot = [];

%Simlulation
for i = 1:10                 %max number of steps allowed (depends on tfinal)
    %Solve until the first terminal event.
    [t,x] = ode45(@(t,x) dx_vector_field(t,x,a), [tstart tfinal], x_plus, options);
    nt = length(t);
    
    %Setting the new initial conditions based on impact map
    x_plus = impact_map(x(nt,:));
    
    x_plus = x_plus(1:6);       % Only positions and velocities needed as initial conditions
    
    %
    %Stop conditions
    if tstart >= tfinal || nt-refine < 1
        break
    end
    %}
    
    % A good guess of a valid first timestep is the length of the last valid
    % timestep, so use it for faster computation.  'refine' is 4 by default.
    options = odeset(options,'InitialStep',t(nt)-t(nt-refine),'MaxStep',t(nt)-t(1));
    tstart = t(nt);
    
    %Save data
    time_out = [time_out; t(2:nt)];
    states = [states;x(2:nt,:)];%size([dhdq; D(1,1:3)])
    foot = [foot; (x(1,1)+x(end,1))/2*ones(length(x(:,1)),1) + states(end,1)  - x(end,1)]; %remember foot location
    
end

%
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

%Event function
    function [limits,isterminal,direction] = events(~,z)
        
        q1 = z(1);
        s = (z_plus(1) - q1)/delq;   %normalized general coordinate
        
        limits(1) = s-1; 	%check when gait is at the end
        limits(2) = s+0.1;      %check if gait rolls back to 0
        isterminal = [1 1];    	% Halt integation
        direction = [];       	%The zero can be approached from either direction
    end

    function dx = dx_vector_field(~,x,a)
        
        %Computes vector field x_dot = f(x) + g(x)*u
        
        %Required inputs: x - all states [q; q_dot] and a - Bezier coeffs
        
        t_temp = num2cell(x(1:3)); [q1, q2, q3] = t_temp{:};
        dt_temp = num2cell(x(4:6)); [dq1, dq2, dq3] = dt_temp{:};
        
        [D,C,G,B] = state_matrix(x);        %Get state matrix using full states x
        
        %Compute fx and gx
        Fx = [x(4:6); D\(-C*x(4:6)-G)];
        Gx = [zeros(3,2);D\B];
        
        % Bezier coefficients for q2
        a21 = a(1); a22 = a(2); a23 = a(3); a24 = a(4); a25 = a(5);
        a2 = a(1:5);
        % Bezier coefficients for q3
        a31 = a(6); a32 = a(7); a33 = a(8); a34 = a(9); a35 = a(10);
        a3 = a(6:10);
        
        %s variable used in bezier polynomial
        %delq = deg2rad(30);             %difference between min and max q1
        s = (z_plus(1) - q1)/delq;         %normalized general coordinate
        
        b2 = bezier(s,4,a2);
        b3 = bezier(s,4,a3);
        
        %Output variable
        h = [q2 - b2; q3 - b3];         %y = h(x) = Hq - hd
        
        dh_dx = [ 0, 1, 0, 0, 0, 0;...
            0, 0, 1, 0, 0, 0];
        dh_dx(1,1) = ((10000*a22*((2500*q1)/1309 - 3/2)^3)/1309 - (10000*a21*((2500*q1)/1309 - 3/2)^3)/1309 + (10000*a24*((2500*q1)/1309 - 1/2)^3)/1309 - (10000*a25*((2500*q1)/1309 - 1/2)^3)/1309 + (30000*a22*((2500*q1)/1309 - 1/2)*((2500*q1)/1309 - 3/2)^2)/1309 + (30000*a24*((2500*q1)/1309 - 1/2)^2*((2500*q1)/1309 - 3/2))/1309 - 6*a23*((2500*q1)/1309 - 3/2)^2*((12500000*q1)/1713481 - 2500/1309) - 6*a23*((2500*q1)/1309 - 1/2)^2*((12500000*q1)/1713481 - 7500/1309))/delq;
        dh_dx(2,1) = ((10000*a32*((2500*q1)/1309 - 3/2)^3)/1309 - (10000*a31*((2500*q1)/1309 - 3/2)^3)/1309 + (10000*a34*((2500*q1)/1309 - 1/2)^3)/1309 - (10000*a35*((2500*q1)/1309 - 1/2)^3)/1309 + (30000*a32*((2500*q1)/1309 - 1/2)*((2500*q1)/1309 - 3/2)^2)/1309 + (30000*a34*((2500*q1)/1309 - 1/2)^2*((2500*q1)/1309 - 3/2))/1309 - 6*a33*((2500*q1)/1309 - 3/2)^2*((12500000*q1)/1713481 - 2500/1309) - 6*a33*((2500*q1)/1309 - 1/2)^2*((12500000*q1)/1713481 - 7500/1309))/delq;
        
        Lfh = dh_dx*Fx;
        
        dLfh = [ -dq1*((75000000*a21*((2500*q1)/1309 - 3/2)^2)/1713481 + (75000000*a23*((2500*q1)/1309 - 1/2)^2)/1713481 - (150000000*a22*((2500*q1)/1309 - 3/2)^2)/1713481 - (150000000*a24*((2500*q1)/1309 - 1/2)^2)/1713481 + (75000000*a23*((2500*q1)/1309 - 3/2)^2)/1713481 + (75000000*a25*((2500*q1)/1309 - 1/2)^2)/1713481 - (30000*a24*((2500*q1)/1309 - 3/2)*((12500000*q1)/1713481 - 2500/1309))/1309 - (30000*a22*((2500*q1)/1309 - 1/2)*((12500000*q1)/1713481 - 7500/1309))/1309 + 12*a23*((12500000*q1)/1713481 - 2500/1309)*((12500000*q1)/1713481 - 7500/1309)), 0, 0, (10000*a22*((2500*q1)/1309 - 3/2)^3)/1309 - (10000*a21*((2500*q1)/1309 - 3/2)^3)/1309 + (10000*a24*((2500*q1)/1309 - 1/2)^3)/1309 - (10000*a25*((2500*q1)/1309 - 1/2)^3)/1309 + (30000*a22*((2500*q1)/1309 - 1/2)*((2500*q1)/1309 - 3/2)^2)/1309 + (30000*a24*((2500*q1)/1309 - 1/2)^2*((2500*q1)/1309 - 3/2))/1309 - 6*a23*((2500*q1)/1309 - 3/2)^2*((12500000*q1)/1713481 - 2500/1309) - 6*a23*((2500*q1)/1309 - 1/2)^2*((12500000*q1)/1713481 - 7500/1309), 1, 0;...
            -dq1*((75000000*a31*((2500*q1)/1309 - 3/2)^2)/1713481 + (75000000*a33*((2500*q1)/1309 - 1/2)^2)/1713481 - (150000000*a32*((2500*q1)/1309 - 3/2)^2)/1713481 - (150000000*a34*((2500*q1)/1309 - 1/2)^2)/1713481 + (75000000*a33*((2500*q1)/1309 - 3/2)^2)/1713481 + (75000000*a35*((2500*q1)/1309 - 1/2)^2)/1713481 - (30000*a34*((2500*q1)/1309 - 3/2)*((12500000*q1)/1713481 - 2500/1309))/1309 - (30000*a32*((2500*q1)/1309 - 1/2)*((12500000*q1)/1713481 - 7500/1309))/1309 + 12*a33*((12500000*q1)/1713481 - 2500/1309)*((12500000*q1)/1713481 - 7500/1309)), 0, 0, (10000*a32*((2500*q1)/1309 - 3/2)^3)/1309 - (10000*a31*((2500*q1)/1309 - 3/2)^3)/1309 + (10000*a34*((2500*q1)/1309 - 1/2)^3)/1309 - (10000*a35*((2500*q1)/1309 - 1/2)^3)/1309 + (30000*a32*((2500*q1)/1309 - 1/2)*((2500*q1)/1309 - 3/2)^2)/1309 + (30000*a34*((2500*q1)/1309 - 1/2)^2*((2500*q1)/1309 - 3/2))/1309 - 6*a33*((2500*q1)/1309 - 3/2)^2*((12500000*q1)/1713481 - 2500/1309) - 6*a33*((2500*q1)/1309 - 1/2)^2*((12500000*q1)/1713481 - 7500/1309), 0, 1];
        
        epsilon = 0.1; alp = 0.9;
        
        %scaling
        Lfh = epsilon*Lfh;
        
        phi(1) = h(1) + 1/(2 - alp)*sign(Lfh(1))*abs(Lfh(1))^(2-alp);
        phi(2) = h(2) + 1/(2 - alp)*sign(Lfh(2))*abs(Lfh(2))^(2-alp);
        
        psi(1,1) = -sign(Lfh(1))*abs(Lfh(1))^alp - sign(phi(1))*abs(phi(1))^(alp/(2-alp));
        psi(2,1) = -sign(Lfh(2))*abs(Lfh(2))^alp - sign(phi(2))*abs(phi(2))^(alp/(2-alp));
        
        v = 1/epsilon^2*psi;
        
        %{
%PD control
q2_d = deg2rad(30); q3_d = deg2rad(180); dq2_d = 2; dq3_d = 1;
e = [q2 - q2_d; q3 - q3_d]; de = [dq2 - dq2_d; dq3 - dq3_d];
Kp = [1000, 0; 0 1000];
Kd = [1000, 0; 0 1000];
v = -Kp*e - Kd*de;

%v = -[2000; 2000];
        %}
        
        %u = LgLfh^-1*(v - L2fh) needed to cancel non-linearity
        u = (dLfh*Gx)\(v - dLfh*Fx);
        
        dx = Fx + Gx*u;
        
    end

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