% optimize leg params
%
% by Alireza Ramezani, 5-21-2018, Pasadena, CA

close all; clear all; clc;

% leg params
func_params;


% optimize leg params
[x_sol,fval,exitflag,output,lambda,grad,hessian] = func_opt_leg_dims(q1, q2, R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4, x0(1), x0(2), x0(3), x0(4), x0(5), x0(6), fy);