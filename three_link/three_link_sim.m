
global stride_angle
stride_angle = pi/9;                        %Stride angle 
x0 = [-stride_angle; pi/6; pi/6; 2; 2; 2];       %I.C. [thetas; velocities]
tstart = 0; tfinal = 10;                    %max time per swing

refine = 4;
options = odeset('Events',@events,'OutputFcn',@odeplot,'OutputSel',1,'Refine',refine);

fig = figure;
vel = []; x_pos = []; time = []; dist = []; sum = []; foot = [];

for i = 1:5                 %number of steps
    % Solve until the first terminal event.
    [t,x] = ode45(@(t,x) state_eqn_3link(t,x), [tstart tfinal], x0, options);
   if ~ishold
      hold on
   end
      
   % Setting the new initial conditions based on impact map
   x0 = impact_model_3link(x);
      
   nt = length(t);
   
   % A good guess of a valid first timestep is the length of the last valid
   % timestep, so use it for faster computation.  'refine' is 4 by default.
   options = odeset(options,'InitialStep',t(nt)-t(nt-refine),'MaxStep',t(nt)-t(1));
   tstart = t(nt);
   
   vel = [vel; x(:,2)]; x_pos = [x_pos; x(:,1)];
   time = [time; t];
   
   %modifying position to save
   if i == 1 
       dist = x(:,1);
   else
       sum = x(:,1) + abs(x(1,1)) + dist(end,1);
       dist = [dist; sum];
   end
   foot = [foot; (x(1,1)+x(end,1))/2*ones(length(x(:,1)),1) + dist(end,1)  - x(end,1)]; %remember foot location
   
   if x0(2) <= 0
       break                    %stop if velocity is <= 0
   end
end

function [limits,isterminal,direction] = events(t,x)

global stride_angle

limits = stride_angle - x(1);       %check when stance leg reaches desired angle
isterminal = 1;                     % Halt integration
direction = [-1,-1];                % The zero can be approached from either direction

end