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

% optimized values from 2018_5_22_2_46_32.596.mat
q1 = 1.7453;
q2 = 2.7925;
R1 = 0.0711;
R2 = 0.3492;
R3 = 0.1347;
Rc1 = 0.0711;
Rc2 = 0.2080;
Rc3 = 0.2602;
e1 = 0.0700;
e2 = 0.1000;
e3 = 0.0300;
e4 = 0.0500;

% these are the exact solutions from the solver for the above dims
x0(1) = 108.924;
x0(2) = 119.1635;
x0(3) = 31.9093;
x0(4) = 59.9225;
x0(5) = 123.8391;
x0(6) = 57.0749;
x0 = x0.*pi/180;

% load optimized params
% load('optsols\kin\2018_5_22_1_50_41.1.mat');
% load('optsols\kin\2018_5_22_2_33_12.727.mat');
% load('optsols\kin\2018_5_22_2_46_32.596.mat');


% force "F" acts at point p9
fx = 0; % x component [N]

% Note: this is the approximation of the total weight I made on 5-16-2018
fy = 30; % y component [N] 

F = [fx, fy].';