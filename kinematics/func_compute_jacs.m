% compute jacobian terms
% Note: please see "generate_constraints.m" (or my notes).
% By Alireza Ramezani, 5-16-2018, Pasadena, CA
% 
% x:
% alpha1
% alpha2
% alpha3
% alpha4
% alpha5
% alpha6

% motor angles:
% q1
% q2

% params:
% R1 R2 R3 Rc1 Rc2 Rc3 e1 e2 e3 e4

% output:
% J1: jacobian(p9,[q1;q2])
% J2: jacobian(p9,[alpha1;alpha2;alpha3;alpha4;alpha5;alpha6])
% Jc1: jacobian("holonomic constraints vector",[q1;q2])
% Jc2: jacobian("holonomic constraints vector",[alpha1;alpha2;alpha3;alpha4;alpha5;alpha6])

function [J1,J2,Jc1,Jc2] = func_compute_jacs(x, q1, q2, R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4)

% alphas...
alpha1 = x(1);
alpha2 = x(2);
alpha3 = x(3);
alpha4 = x(4);
alpha5 = x(5);
alpha6 = x(6);

% J1:
J1 = [ 0, e2*sin(alpha1 + alpha2 - q2) - R2*sin(alpha1 + alpha2 - q2) - R3*sin(alpha1 + alpha2 + alpha4 - q2) - Rc1*sin(q2) + Rc2*sin(alpha1 - q2) - e3*sin(alpha1 - q2);...
0, R2*cos(alpha1 + alpha2 - q2) - e2*cos(alpha1 + alpha2 - q2) + R3*cos(alpha1 + alpha2 + alpha4 - q2) - Rc1*cos(q2) - Rc2*cos(alpha1 - q2) + e3*cos(alpha1 - q2)];

% J2:
J2 = [R2*sin(alpha1 + alpha2 - q2) - e2*sin(alpha1 + alpha2 - q2) + R3*sin(alpha1 + alpha2 + alpha4 - q2) - Rc2*sin(alpha1 - q2) + e3*sin(alpha1 - q2), R2*sin(alpha1 + alpha2 - q2) - e2*sin(alpha1 + alpha2 - q2) + R3*sin(alpha1 + alpha2 + alpha4 - q2), 0,  R3*sin(alpha1 + alpha2 + alpha4 - q2), 0, 0;...
e2*cos(alpha1 + alpha2 - q2) - R2*cos(alpha1 + alpha2 - q2) - R3*cos(alpha1 + alpha2 + alpha4 - q2) + Rc2*cos(alpha1 - q2) - e3*cos(alpha1 - q2), e2*cos(alpha1 + alpha2 - q2) - R2*cos(alpha1 + alpha2 - q2) - R3*cos(alpha1 + alpha2 + alpha4 - q2), 0, -R3*cos(alpha1 + alpha2 + alpha4 - q2), 0, 0];

% Jc1 (holonomic constraints):
Jc1 = [  0, e2*sin(alpha1 + alpha2 - q2) + R1*sin(alpha1 + alpha2 + alpha3 - q2) - Rc1*sin(q2) + Rc2*sin(alpha1 - q2) - e3*sin(alpha1 - q2);...
0, e3*cos(alpha1 - q2) - R1*cos(alpha1 + alpha2 + alpha3 - q2) - Rc1*cos(q2) - Rc2*cos(alpha1 - q2) - e2*cos(alpha1 + alpha2 - q2);...
-1,                                                                                                                              -1;...
0,                                                                                                                               0;...
0,                                                                                                                               0;...
0,                                                                                                                               0];

% Jc2 (holonomic constraints):
Jc2 = [e3*sin(alpha1 - q2) - R1*sin(alpha1 + alpha2 + alpha3 - q2) - Rc2*sin(alpha1 - q2) - e2*sin(alpha1 + alpha2 - q2), - e2*sin(alpha1 + alpha2 - q2) - R1*sin(alpha1 + alpha2 + alpha3 - q2), -R1*sin(alpha1 + alpha2 + alpha3 - q2),                                                                              0,                                                             0,                                 0;...
e2*cos(alpha1 + alpha2 - q2) + R1*cos(alpha1 + alpha2 + alpha3 - q2) + Rc2*cos(alpha1 - q2) - e3*cos(alpha1 - q2),   e2*cos(alpha1 + alpha2 - q2) + R1*cos(alpha1 + alpha2 + alpha3 - q2),  R1*cos(alpha1 + alpha2 + alpha3 - q2),                                                                              0,                                                             0,                                 0;...
                                                                                                                 1,                                                                      1,                                      1,                                                                              0,                                                             0,                                 0;...
                                                                                                                 0,                                                                      0,                                      0,   Rc3*sin(alpha4 + alpha5) + e4*sin(alpha4) + e3*sin(alpha4 + alpha5 + alpha6),   Rc3*sin(alpha4 + alpha5) + e3*sin(alpha4 + alpha5 + alpha6),  e3*sin(alpha4 + alpha5 + alpha6);...
                                                                                                                 0,                                                                      0,                                      0, - Rc3*cos(alpha4 + alpha5) - e4*cos(alpha4) - e3*cos(alpha4 + alpha5 + alpha6), - Rc3*cos(alpha4 + alpha5) - e3*cos(alpha4 + alpha5 + alpha6), -e3*cos(alpha4 + alpha5 + alpha6);...
                                                                                                                 0,                                                                      1,                                      0,                                                                              1,                                                             1,                                 1];

end
