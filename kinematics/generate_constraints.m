clear all; close all;

% define parameters
syms R1 R2 R3 Rc1 Rc2 Rc3 e1 e2 e3 e4
syms q1 q2 alpha1 alpha2 alpha3 alpha4 alpha alpha5 alpha6
syms p0 p1 p2 p3 p4 p5 p6 p7 p8 p9
syms v1 v2 v3 v4 v5 v6 v7 v8 v9

%  
q = [q1,q2].';
q_pr = [alpha1,alpha2,alpha3,alpha4,alpha5,alpha6].';

% define transformation matrix
T = inline('[1 0 0; x cos(q) -sin(q);0 sin(q) cos(q)]','x','q');

disp('generating constraints ...');

% first loop in the leg mechanism
T10 = T(-e1,2*pi-q2);
T21 = T(Rc1,alpha1);
T32 = T(Rc2-e3,alpha2);
T43 = T(e2,alpha3);
T04 = T(R1,-q1);
tmp = T10*T21*T32*T43*T04;
% tmp = T10*T21*T32-eye(3)/(T43*T04);
T1 = simplify(tmp);

% please look at T matrix to understand why I picked this const...
f1 = [T1(2,1);... % x-position
        T1(3,1);...  % y-position
            alpha1 + alpha2 + alpha3 - q1 - q2];  % sum of angles in polygon

% second loop in the leg mechanism
T53 = T(-(R2-e2),alpha4+pi);
T65 = T(e4,alpha5);
T95 = T(R3,0);
T76 = T(Rc3,alpha6);
T37 = T(e3,alpha2-pi);
tmp = T53*T65*T76*T37;
T2 = simplify(tmp);

% please look at T matrix to understand why I picked this const...
f2 = [T2(2,1);... % x-position
        T2(3,1);...  % y-position
            alpha2 + alpha4 + alpha5 + alpha6 - 2*pi];  % sum of angles in polygon

        
disp('generating points ...');
p1 = [1,0,0].';

tmp = T10;
p2 = simplify(tmp(2:3,1));

tmp = T10*T21;
p3 = simplify(tmp(2:3,1));

tmp = T10*T21*T32;
p4 = simplify(tmp(2:3,1));

tmp = T10*T21*T32*T43;
p5 = simplify(tmp(2:3,1));

tmp = T10*T21*T32*T53;
p6 = simplify(tmp(2:3,1));

tmp = T10*T21*T32*T53*T65;
p7 = simplify(tmp(2:3,1));

tmp = T10*T21*T32*T53*T65*T76;
p8 = simplify(tmp(2:3,1));

tmp = T10*T21*T32*T53*T95;
p9 = simplify(tmp(2:3,1));


% v1 = [0,0].';
% v2 = jacobian(p2,q)*dq;
% v3 = jacobian(p3,q)*dq;
% v4 = jacobian(p4,q)*dq;
% v5 = jacobian(p5,q)*dq;
% v6 = jacobian(p6,q)*dq;
% v7 = jacobian(p7,q)*dq;
% v8 = jacobian(p8,q)*dq;
% v9 = jacobian(p9,q)*dq;
% v = [v1;v2;v3;v4;v5;v6;v7;v8;v9];

% Note: I will be using these jacobians to compute the joint torques in the
% closed kinematic chain subject to external ground reaction force acting
% at point p9.
disp('generating jacobians...');
J1 = jacobian(p9,q);
J2 = jacobian(p9,q_pr);

% stacking holonomic constraints
c = [f1;f2];
Jc1 = jacobian(c,q);
Jc2 = jacobian(c,q_pr);



































