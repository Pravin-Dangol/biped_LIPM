% returns fwd kinematics terms (pos and velocities of joints in LEO's leg
% mechanism, p1, p2, etc.)
% 12/22/2017 Writen by Alirza Ramezani, Pasadena, CA
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
%
% output:
% p1, p2 ,p3, ... p9
% v1, v2 ,v3, ... v9
function [p1,p2,p3,p4,p5,p6,p7,p8,p9,...
    v1,v2,v3,v4,v5,v6,v7,v8,v9] = func_leg_fwd_kin(x, q1, q2, R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4)

alpha1 = x(1);
alpha2 = x(2);
alpha3 = x(3);
alpha4 = x(4);
alpha5 = x(5);
alpha6 = x(6);

p1 = [0,0].';   %p1

p2 = [-e1,0].';

p3 = [Rc1*cos(q2) - e1, -Rc1*sin(q2)].';

p4 = [Rc1*cos(q2) - e1 + cos(alpha1 - q2)*(Rc2 - e3),...
      sin(alpha1 - q2)*(Rc2 - e3) - Rc1*sin(q2)].';

p5 = [e2*cos(alpha1 + alpha2 - q2) - e1 + Rc1*cos(q2) + Rc2*cos(alpha1 - q2) - e3*cos(alpha1 - q2),...
      e2*sin(alpha1 + alpha2 - q2) - Rc1*sin(q2) + Rc2*sin(alpha1 - q2) - e3*sin(alpha1 - q2)].';
  
p6 = [Rc1*cos(q2) - cos(alpha1 + alpha2 - q2)*(R2 - e2) - e1 + cos(alpha1 - q2)*(Rc2 - e3),...
      sin(alpha1 - q2)*(Rc2 - e3) - Rc1*sin(q2) - sin(alpha1 + alpha2 - q2)*(R2 - e2)].';  
  
p7 = [e2*cos(alpha1 + alpha2 - q2) - R2*cos(alpha1 + alpha2 - q2) - e1 - e4*cos(alpha1 + alpha2 + alpha4 - q2) + Rc1*cos(q2) + Rc2*cos(alpha1 - q2) - e3*cos(alpha1 - q2),...
      e2*sin(alpha1 + alpha2 - q2) - R2*sin(alpha1 + alpha2 - q2) - e4*sin(alpha1 + alpha2 + alpha4 - q2) - Rc1*sin(q2) + Rc2*sin(alpha1 - q2) - e3*sin(alpha1 - q2)].';  

p8 = [ e2*cos(alpha1 + alpha2 - q2) - R2*cos(alpha1 + alpha2 - q2) - e1 - e4*cos(alpha1 + alpha2 + alpha4 - q2) + Rc1*cos(q2) + Rc2*cos(alpha1 - q2) - e3*cos(alpha1 - q2) - Rc3*cos(alpha1 + alpha2 + alpha4 + alpha5 - q2),...
      e2*sin(alpha1 + alpha2 - q2) - R2*sin(alpha1 + alpha2 - q2) - e4*sin(alpha1 + alpha2 + alpha4 - q2) - Rc1*sin(q2) + Rc2*sin(alpha1 - q2) - e3*sin(alpha1 - q2) - Rc3*sin(alpha1 + alpha2 + alpha4 + alpha5 - q2)].';

p9 = [e2*cos(alpha1 + alpha2 - q2) - R2*cos(alpha1 + alpha2 - q2) - e1 - R3*cos(alpha1 + alpha2 + alpha4 - q2) + Rc1*cos(q2) + Rc2*cos(alpha1 - q2) - e3*cos(alpha1 - q2),...
      e2*sin(alpha1 + alpha2 - q2) - R2*sin(alpha1 + alpha2 - q2) - R3*sin(alpha1 + alpha2 + alpha4 - q2) - Rc1*sin(q2) + Rc2*sin(alpha1 - q2) - e3*sin(alpha1 - q2)].';


v1 = [];
v2 = [];
v3 = [];
v4 = [];
v5 = [];
v6 = [];
v7 = [];
v8 = [];
v9 = [];



end