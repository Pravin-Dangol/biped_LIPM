% map small variations in leg-end to small variations in joint angles at
% points q = (q1;q2) and qp = (qp1;qp2;qp3;qp4;qp5;qp6)
%
% input:
%   q(1): HD1 angle
%   q(2): HD2 angle
%
%   qp(1): alpha1 (look at notes)
%   qp(2): alpha2
%   qp(3): alpha3
%   qp(4): alpha4
%   qp(5): alpha5
%   qp(6): alpha6
%
%   delta_p(1): horizontal variations in leg-end position  
%   delta_p(2): vertical variations in leg-end position
%
%   params: R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4
%
% output:
%   delta_q(1): variations in HD1 angle
%   delta_q(2): variations in HD1 angle
%
% by Alireza Ramezani, 10-6-2018, Pawtucket, RI
function [delta_q] = func_task_space_to_joint_space(q,qp,delta_p,R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4)

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
[J1,J2,Jc1,Jc2] = func_compute_jacs(qp, q(1), q(2), R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4);

J = [J1,J2;Jc1,Jc2];

tmp = J\[delta_p;zeros(6,1)];

delta_q = tmp(1:2);
end