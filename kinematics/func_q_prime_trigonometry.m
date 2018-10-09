% solve q primes using simple trigonometric rules
%
% by Alireza Ramezani, 10-3-2018, pawtucket, RI
%
% Inputs:
%   q1
%   q2
%
% params:
%   R1 R2 R3 Rc1 Rc2 Rc3 e1 e2 e3 e4
%
% Output: (see my notes to see the def. of these angles)
%   qp1
%   qp2
%   qp3
%   qp4
%   qp5
%   qp6
function [qp1,qp2,qp3,qp4,qp5,qp6] = func_q_prime_trigonometry(q1, q2, R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4)

l1 = sqrt(e1^2+Rc1^2-2*e1*Rc1*cos(q2));

% 3rd output is q2
[alpha1,alpha2,~] = func_triangle_angles(l1,e1,Rc1);

l2 = sqrt(l1^2+R1^2-2*l1*R1*cos(q1-alpha1));

% 3rd output is q1-alpha2
[alpha3,alpha4,~] = func_triangle_angles(l2,R1,l1);

% 2nd output is pi-qp2
[alpha5,tmp,alpha6] = func_triangle_angles(Rc2-e3,e2,l2);

l3 = sqrt((R2-e2)^2+e3^2-2*(R2-e2)*e3*cos(tmp));

% 3rd output is pi-qp2
[alpha7,alpha8,~] = func_triangle_angles(l3,e3,R2-e2);

[alpha9,alpha11,alpha10] = func_triangle_angles(e4,Rc3,l3);


qp1 = pi-(alpha1+alpha3+alpha5);
qp2 = pi-tmp;
qp3 = pi-(alpha4+alpha6);
qp4 = pi-(alpha7+alpha9);
qp5 = pi-alpha11;
qp6 = pi-(alpha8+alpha10);

% use trigonometry to obtain the angles of a triangle with
%
% input:
%   lines with dims a,b and c
%
% output:
%   theta1: angle between a and c
%   theta2: angle between c and b
%   theta3: angle between b and a
    function [theta1,theta2,theta3] = func_triangle_angles(a,b,c)
        c1 = func_cos_theta(b,c,a);
        c2 = func_cos_theta(c,a,b);
        c3 = func_cos_theta(a,b,c);
        
        theta1 = acos(c1);
        theta2 = acos(c2);
        theta3 = acos(c3);
    end

% obtain cos(theta) based on given params a, b, c dims
%
% input:
%   a: dim of line in front of theta
%   b, c: dim of lines adjacent to theta
%
% outout:
%   cos(theta)
    function q = func_cos_theta(a,b,c)
        q = (a^2-b^2-c^2)/(-2*b*c);
    end

end