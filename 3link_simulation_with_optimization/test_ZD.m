%Simulation of a single step of ZD to compare with full dynamics

f = [-0.2018    0.0463    0.4108   -0.0212    1.7114    0.5025   -1.9435    0.4918  0.5809];
M = 4;
alpha2 = [f(5), 2*f(5)-f(4), f(3:5)];
alpha3 = [f(8), 2*f(8)-f(7), f(6:8)];

global a
a = [alpha2, alpha3];

x0 = map_z_to_x([f(1), f(2)],a);   %I.C. [thetas; velocities]

x_plus = impact_map(x0);
x_plus = x_plus(1:6)';

z_plus = [x_plus(1), x_plus(4)];

tstart = 0; tfinal = 50;                    %max time per swing
time_out = tstart;

refine = 4; options = odeset('Events',@events,'Refine',refine,'RelTol',10^-5,'AbsTol',10^-6);

[t,z] = ode45(@(t,z) ZD_states(t,z,a), [tstart tfinal], z_plus, options);

nt = length(t); time_out = [time_out; t(2:nt)];

[m,~] = size(z); states = zeros(m,6);
for i = 1:m
    states(i,:) = map_z_to_x(z(i,:),a);
end

%
%Plots
figure(2)
subplot(2,1,1)
plot(time_out,states(:,1),time_out,states(:,2),'-.',time_out,states(:,3),'--')
legend('\theta_1','\theta_2','\theta_3')
title('Joint angles')
subplot(2,1,2)
plot(time_out,states(:,4),time_out,states(:,5),'-.',time_out,states(:,6),'--')
legend('$\dot{\theta_1}$','$\dot{\theta_2}$','$\dot{\theta_3}$','Interpreter','latex')
title('Joint velocities')
%}

%Event function
function [limits,isterminal,direction] = events(~,z)
global a
[r,~,~,~,~,~] = model_params_3link;
q1d = control_params_3link;

q1 = z(1);
x = map_z_to_x(z,a);
q2 = x(2);

limits(1) = r*cos(q1) + r*cos(q1 - q2); 	%check when stance leg reaches desired angle
%limits(2) = r*sin(z(1)) - 0.1*r;    %check if leg is close to ground
isterminal = 1;                     % Halt integation
direction = [];                      %The zero can be approached from either direction

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

delq = deg2rad(30);              %difference between min and max q1
s = (q1 + delq/2)/delq;   %normalized general coordinate

dLsb2 = -dq1/delq*(3*s^2*(4*a24 - 4*a25) - s^2*(12*a23 - 12*a24) - 3*(s - 1)^2*(4*a21 - 4*a22) +...
    (s - 1)^2*(12*a22 - 12*a23) - 2*s*(s - 1)*(12*a23 - 12*a24) + s*(2*s - 2)*(12*a22 - 12*a23));

dLsb3 = -dq1/delq*(3*s^2*(4*a34 - 4*a35) - s^2*(12*a33 - 12*a34) - 3*(s - 1)^2*(4*a31 - 4*a32) +...
    (s - 1)^2*(12*a32 - 12*a33) - 2*s*(s - 1)*(12*a33 - 12*a34) + s*(2*s - 2)*(12*a32 - 12*a33));

beta1 = [dLsb2; dLsb3]*dq1/delq;

%beta1 = [(6250000*dq1^2*((7500*((2500*q1)/1309 - 3/2)^2*(4*a21 - 4*a22))/1309 - (7500*((2500*q1)/1309 - 1/2)^2*(4*a24 - 4*a25))/1309 + (2500*((2500*q1)/1309 - 1/2)^2*(12*a23 - 12*a24))/1309 - (2500*((2500*q1)/1309 - 3/2)^2*(12*a22 - 12*a23))/1309 + ((2500*q1)/1309 - 3/2)*((12500000*q1)/1713481 - 2500/1309)*(12*a23 - 12*a24) - ((2500*q1)/1309 - 1/2)*((12500000*q1)/1713481 - 7500/1309)*(12*a22 - 12*a23)))/1713481;...
%    (6250000*dq1^2*((7500*((2500*q1)/1309 - 3/2)^2*(4*a21 - 4*a22))/1309 - (7500*((2500*q1)/1309 - 1/2)^2*(4*a24 - 4*a25))/1309 + (2500*((2500*q1)/1309 - 1/2)^2*(12*a23 - 12*a24))/1309 - (2500*((2500*q1)/1309 - 3/2)^2*(12*a22 - 12*a23))/1309 + ((2500*q1)/1309 - 3/2)*((12500000*q1)/1713481 - 2500/1309)*(12*a23 - 12*a24) - ((2500*q1)/1309 - 1/2)*((12500000*q1)/1713481 - 7500/1309)*(12*a22 - 12*a23)))/1713481];

db_ds2 = d_ds_bezier(s,4,a_2); db_ds3 = d_ds_bezier(s,4,a_3);

beta2 = [db_ds2; db_ds3]/delq;

dz(1) = z(2);
dz(2) = (D1 + D2*beta2)\(-D2*beta1 - H1);
%dz(2) = -G(1,1);

dz = [dz(1), dz(2)]';

end