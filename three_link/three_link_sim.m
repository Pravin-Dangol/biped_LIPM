
global stride_angle
stride_angle = pi/8;                        %Stride angle
x0 = [-stride_angle; stride_angle; pi/9; 1; 1; 1];       %I.C. [thetas; velocities]
tstart = 0; tfinal = 10;                    %max time per swing

refine = 4; options = odeset('Events',@events,'OutputSel',1,'Refine',refine);    %'OutputFcn',@odeplot,

time_out = tstart; states = x0.'; foot = [];

%Simlulation
for i = 1:1000                 %number of steps
    %Solve until the first terminal event.
    [t,x] = ode45(@(t,x) state_eqn_3link(t,x), [tstart tfinal], x0, options);
    nt = length(t);

    %Setting the new initial conditions based on impact map
    x0 = impact_map_3link(x(nt,:));
    
    % Only positions and velocities needed as initial conditions
    x0 = x0(1:6);
    
    % A good guess of a valid first timestep is the length of the last valid
    % timestep, so use it for faster computation.  'refine' is 4 by default.
    options = odeset(options,'InitialStep',t(nt)-t(nt-refine),'MaxStep',t(nt)-t(1));
    tstart = t(nt);
    
    %Save data
    time_out = [time_out; t(2:nt)]; states = [states;x(2:nt,:)];
    foot = [foot; (x(1,1)+x(end,1))/2*ones(length(x(:,1)),1) + states(end,1)  - x(end,1)]; %remember foot location
    
    %Stop copnditions
    if x0(4) <= 0 
        disp(['Velocity drop at step number: ', num2str(i)-1])
        break                	%stop if velocity is <= 0
    elseif cos(x(1)) <= 0.1 || cos(x(2)) <= 0.1 
        disp(['about to fall at step number: ', num2str(i)-1])
        break                   %stop if leg is too close to ground
    elseif tstart >= tfinal
        break
    end
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

function dx = state_eqn_3link(t,x)

t_temp = num2cell(x(1:3)); [t1, t2, t3] = t_temp{:};

[D,C,G,B] = state_matrix_3link(t,x);

fx = D\(-C*x(4:6)-G);   gx = D\B;

t3d = pi/12;         %arbitrary desired torso angle
t1d = pi/8;
%set torso to desired angle, assure two legs are symmetric (t1 = -t2)
y = [t1 - t1d + t2; t3 - t3d];    %y = h(x)
dh_dx = [1, 1, 0;...
    0, 0, 1];       %dy/dx

K = [2, 0;...
    -1.4, 2.5];
%[2, 0; -1.4, 2.5]; %stumbled for 567 steps
%[5, 0; 0, 12]; %stumbled for 42 steps

u = (dh_dx*gx)\(-K*y - dh_dx*fx); %needed to cancel non-linearity

dx(1:3) = x(4:6);   dx(4:6) = fx + gx*u;
dx = dx';

end

%Event function
function [limits,isterminal,direction] = events(~,x)

global stride_angle
[r,~,~,~,~,~] = model_params_3link;

limits(1) = stride_angle - x(1);       %check when stance leg reaches desired angle
limits(2) = r*cos(x(1)) - 0.1*r;    %check if leg is close to ground
isterminal = [1, 1];                     % Halt integation
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