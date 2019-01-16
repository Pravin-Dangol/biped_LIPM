%Simulates full hybrid dynamics with parametrized output function

%Output of fmincon [q1, q1_dot, alpha,J]
f = [0.1773    2.4665   -1.5124    1.2887    1.5019    0.4561   -2.9580    1.6317  999.9990];
M = 4; delq = deg2rad(30);              %difference between min and max q1

z_minus = f(1:2); 
%alpha is the Bezier coeffs for q2, gamma for q3
alpha = f(3:5);         %alpha 3 - 5 for q2
gamma = f(6:8);         %alpha 3 - 5 for q3

%Mapping from z to x on boundary pg 140, 141
%At the boundry (right before impact, q1_minus)
q2_minus = alpha(3);            %q- = alpha(0) = alpha(M) (desired)
q3_minus = gamma(3);

a0 = alpha(3); g0 = gamma(3);       %alpha(0) = alpha(M) (desired)

%qb_dot = db/ds*s_dot = M*(alpha(4) - alpha(3))*theta_dot-/(delta_theta)
dq2_minus = M*(alpha(3) - alpha(2))*z_minus(2)/delq;
dq3_minus = M*(gamma(3) - gamma(2))*z_minus(2)/delq;
%Enforcing slope for the 1st 2 coefficients to be equal to qb_dot:
%M*(alpha(1) - alpha(0))*theta_dot-/(delta_theta) = qb_dot, so:
a1 = dq2_minus*delq/(M*z_minus(2)) + a0;
g1 = dq3_minus*delq/(M*z_minus(2)) + g0;

a = [a0, a1, alpha, g0, g1, gamma];

%I.C. [thetas; velocities]
x0 = [z_minus(1), q2_minus, q3_minus, z_minus(2), dq2_minus, dq3_minus];

%Apply impact map
x_plus = impact_map(x0);
x_plus = x_plus(1:6)';

tstart = 0; tfinal = 10;                    %max time per swing

refine = 4; options = odeset('Events',@events,'Refine',refine);    %'OutputFcn',@odeplot,'OutputSel',1,

time_out = tstart; states = x_plus.'; foot = [];

%Simlulation
for i = 1:1                 %max number of steps allowed (depends on tfinal)
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
function [limits,isterminal,direction] = events(~,x)

[r,~,~,~,~,~] = model_params_3link;
q1 = x(1); q2 = x(2);

limits = r*cos(q1) + r*cos(q1 - q2);     %check when swing leg is close to ground
%limits(2) = r*sin(q1) - 0.1*r;    %check if hip is close to ground
isterminal = 1;                     % Halt integation
direction = [];                      %The zero can be approached from either direction

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