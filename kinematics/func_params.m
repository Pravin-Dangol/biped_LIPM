% param: from CAD model
% Note: 
%       1) angles [rad] and
%       2) dims [m]
% R1 R2 R3 Rc1 Rc2 Rc3 e1 e2 e3 e4

% motor angles [rad]
% q1 = 140*pi/180;  % q1
% q2 = 110*pi/180;  % q2
% R1 = 0.05;  % R1 [m]
% R2 = 0.25;  % R2 [m]
% R3 = 0.1;  % R3 [m]
% Rc1 = 0.06; % Rc1 [m]
% Rc2 = 0.15; % Rc2 [m]
% Rc3 = 0.12; % Rc3 [m]
% e1 = 0.077;  % e1 [m]
% e2 = 0.13; % e2 [m]
% e3 = 0.05; % e3 [m]
% e4 = 0.040; % e4 [m]

% these are the exact solutions from the solver for the above dims
x0(1) = 56.2106;    % alpha1
x0(2) = 107.6062;   % alpha2
x0(3) = 86.1832;    % alpha3
x0(4) = 68.4871;    % alpha4
x0(5) = 106.5198;   % alpha5
x0(6) = 77.3869;    % alpha6
x0 = x0.*pi/180;

% load optimized params
% load('optsols\kin\2018_5_22_1_50_41.1.mat');
% load('optsols\kin\2018_5_22_2_33_12.727.mat');
load('optsols\kin\2018_5_22_2_46_32.596.mat');


q1 = x_sol(1);
q2 = x_sol(2);
R1 = x_sol(3);
R2 = x_sol(4);
R3 = x_sol(5);
Rc1 = x_sol(6);
Rc2 = x_sol(7);
Rc3 = x_sol(8);
e1 = x_sol(9);
e2 = x_sol(10);
e3 = x_sol(11);
e4 = x_sol(12);



% force "F" acts at point p9
fx = 0; % x component [N]

% Note: this is the approximation of the total weight I made on 5-16-2018
fy = 30; % y component [N] 

F = [fx, fy].';