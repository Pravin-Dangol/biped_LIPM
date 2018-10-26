function dx = state_eqn_3link(t,x)

t_temp = num2cell(x(1:3)); [t1, t2, t3] = t_temp{:};

[D,C,G,B] = state_matrix_3link(t,x);
[r,m,Mh,Mt,l,g0] = model_params_3link;

fx = D\(-C*x(4:6)-G);
gx = D\B;

t3d = pi/12;         %arbitrary desired torso angle
t1d = pi/8;
%set torso to desired angle, assure two legs have opposite angles (t1 = -t2)
y = [t1 - t1d + t2; t3 - t3d];    %y = h(x) 
%dh = [1, 1, 0; 0, 0, 1];       %dy/dx  %jacobian(y,q);
%dLfh = d(dh*fx)/dx;           
dLfh11 = (2*(3*g0*m*cos(t1) + 4*Mt*g0*cos(t1) + 2*Mt*g0*sin(t1 - t3)*sin(t3) + 2*g0*m*sin(t1 - t2)*sin(t2) - 2*Mt*l*t3*t3_d*cos(t1 - t3) + m*r*t2*t2_d*cos(t1 - t3) - 2*Mt*r*t1*t1_d*cos(t1 - t3)^2 + 2*Mt*r*t1*t1_d*sin(t1 - t3)^2 - 2*m*r*t1*t1_d*cos(t1 - t2)^2 + 2*m*r*t1*t1_d*sin(t1 - t2)^2 - 2*Mt*r*t1_d*cos(t1 - t3)*sin(t1 - t3) - 2*m*r*t1_d*cos(t1 - t2)*sin(t1 - t2)))/(r*(4*Mh + 4*Mt + 5*m - 4*Mt*cos(t1 - t3)^2 - 4*m*cos(t1 - t2)^2)) - (2*(8*Mt*g0*sin(t1 - t2)*sin(t1) - 6*g0*m*cos(t1 - t2)*cos(t1) + 6*g0*m*sin(t1 - t2)*sin(t1) + 4*Mh*r*t1_d*sin(t1 - t2) + 4*Mt*r*t1_d*sin(t1 - t2) + 5*m*r*t1_d*sin(t1 - t2) - 8*Mt*g0*cos(t1 - t2)*cos(t1) + 4*Mh*r*t1*t1_d*cos(t1 - t2) + 4*Mt*r*t1*t1_d*cos(t1 - t2) - 4*Mt*r*t1_d*cos(t1 - t3)^2*sin(t1 - t2) - 4*Mt*g0*cos(t1 - t2)*sin(t1 - t3)*sin(t3) - 4*Mt*g0*cos(t1 - t3)*sin(t1 - t2)*sin(t3) + 8*Mt*g0*cos(t1 - t3)*sin(t1 - t3)*sin(t2) + 5*m*r*t1*t1_d*cos(t1 - t2) + 4*Mt*l*t3*t3_d*cos(t1 - t2)*cos(t1 - t3) - 4*Mt*l*t3*t3_d*sin(t1 - t2)*sin(t1 - t3) - 2*m*r*t2*t2_d*cos(t1 - t2)*cos(t1 - t3) + 2*m*r*t2*t2_d*sin(t1 - t2)*sin(t1 - t3) - 4*Mt*r*t1*t1_d*cos(t1 - t2)*sin(t1 - t3)^2 + 4*Mt*r*t1_d*cos(t1 - t2)*cos(t1 - t3)*sin(t1 - t3) + 4*Mt*r*t1*t1_d*cos(t1 - t3)*sin(t1 - t2)*sin(t1 - t3)))/(r*(4*Mh + 4*Mt + 5*m - 4*Mt*cos(t1 - t3)^2 - 4*m*cos(t1 - t2)^2)) + (2*(8*Mt*cos(t1 - t3)*sin(t1 - t3) + 8*m*cos(t1 - t2)*sin(t1 - t2))*(5*g0*m*sin(t2) + 4*Mh*g0*sin(t2) + 4*Mt*g0*sin(t2) - 8*Mt*g0*cos(t1 - t2)*sin(t1) - 6*g0*m*cos(t1 - t2)*sin(t1) - 4*Mt*g0*cos(t1 - t3)^2*sin(t2) + 4*Mt*g0*cos(t1 - t2)*cos(t1 - t3)*sin(t3) + 4*Mh*r*t1*t1_d*sin(t1 - t2) + 4*Mt*r*t1*t1_d*sin(t1 - t2) + 5*m*r*t1*t1_d*sin(t1 - t2) + 4*Mt*l*t3*t3_d*cos(t1 - t2)*sin(t1 - t3) - 2*m*r*t2*t2_d*cos(t1 - t2)*sin(t1 - t3) - 4*Mt*r*t1*t1_d*cos(t1 - t3)^2*sin(t1 - t2) + 4*Mt*r*t1*t1_d*cos(t1 - t2)*cos(t1 - t3)*sin(t1 - t3)))/(r*(- 4*m*cos(t1 - t2)^2 - 4*Mt*cos(t1 - t3)^2 + 4*Mh + 4*Mt + 5*m)^2) + (2*(8*Mt*cos(t1 - t3)*sin(t1 - t3) + 8*m*cos(t1 - t2)*sin(t1 - t2))*(2*Mt*g0*cos(t1 - t3)*sin(t3) - 4*Mt*g0*sin(t1) - 3*g0*m*sin(t1) + 2*g0*m*cos(t1 - t2)*sin(t2) + 2*Mt*l*t3*t3_d*sin(t1 - t3) - m*r*t2*t2_d*sin(t1 - t3) + 2*Mt*r*t1*t1_d*cos(t1 - t3)*sin(t1 - t3) + 2*m*r*t1*t1_d*cos(t1 - t2)*sin(t1 - t2)))/(r*(- 4*m*cos(t1 - t2)^2 - 4*Mt*cos(t1 - t3)^2 + 4*Mh + 4*Mt + 5*m)^2);
dLfh12 = (8*Mt*g0*sin(t1 - t3)*sin(t1) - 6*g0*m*cos(t1 - t3)*cos(t1) + 6*g0*m*sin(t1 - t3)*sin(t1) + 4*Mh*r*t1_d*sin(t1 - t3) + 4*Mt*r*t1_d*sin(t1 - t3) + 5*m*r*t1_d*sin(t1 - t3) - 8*Mt*g0*cos(t1 - t3)*cos(t1) + 4*Mh*r*t1*t1_d*cos(t1 - t3) + 4*Mt*r*t1*t1_d*cos(t1 - t3) + 5*m*r*t1*t1_d*cos(t1 - t3) - 4*m*r*t1_d*cos(t1 - t2)^2*sin(t1 - t3) + 8*g0*m*cos(t1 - t2)*sin(t1 - t2)*sin(t3) - 4*g0*m*cos(t1 - t2)*sin(t1 - t3)*sin(t2) - 4*g0*m*cos(t1 - t3)*sin(t1 - t2)*sin(t2) + 4*Mt*l*t3*t3_d*cos(t1 - t3)^2 - 4*Mt*l*t3*t3_d*sin(t1 - t3)^2 - 2*m*r*t2*t2_d*cos(t1 - t3)^2 + 2*m*r*t2*t2_d*sin(t1 - t3)^2 - 4*m*r*t1*t1_d*cos(t1 - t3)*sin(t1 - t2)^2 + 4*m*r*t1_d*cos(t1 - t2)*cos(t1 - t3)*sin(t1 - t2) + 4*m*r*t1*t1_d*cos(t1 - t2)*sin(t1 - t2)*sin(t1 - t3))/(l*(4*Mh + 4*Mt + 5*m - 4*Mt*cos(t1 - t3)^2 - 4*m*cos(t1 - t2)^2)) - ((8*Mt*cos(t1 - t3)*sin(t1 - t3) + 8*m*cos(t1 - t2)*sin(t1 - t2))*(5*g0*m*sin(t3) + 4*Mh*g0*sin(t3) + 4*Mt*g0*sin(t3) - 8*Mt*g0*cos(t1 - t3)*sin(t1) - 6*g0*m*cos(t1 - t3)*sin(t1) - 4*g0*m*cos(t1 - t2)^2*sin(t3) + 4*Mh*r*t1*t1_d*sin(t1 - t3) + 4*Mt*r*t1*t1_d*sin(t1 - t3) + 4*g0*m*cos(t1 - t2)*cos(t1 - t3)*sin(t2) + 5*m*r*t1*t1_d*sin(t1 - t3) + 4*Mt*l*t3*t3_d*cos(t1 - t3)*sin(t1 - t3) - 2*m*r*t2*t2_d*cos(t1 - t3)*sin(t1 - t3) - 4*m*r*t1*t1_d*cos(t1 - t2)^2*sin(t1 - t3) + 4*m*r*t1*t1_d*cos(t1 - t2)*cos(t1 - t3)*sin(t1 - t2)))/(l*(- 4*m*cos(t1 - t2)^2 - 4*Mt*cos(t1 - t3)^2 + 4*Mh + 4*Mt + 5*m)^2);
dLfh13 = (2*(8*Mt*g0*sin(t1 - t2)*sin(t1) - 4*Mh*g0*cos(t2) - 4*Mt*g0*cos(t2) - 5*g0*m*cos(t2) + 6*g0*m*sin(t1 - t2)*sin(t1) + 4*Mt*g0*cos(t1 - t3)^2*cos(t2) + 4*Mh*r*t1*t1_d*cos(t1 - t2) + 4*Mt*r*t1*t1_d*cos(t1 - t2) - 4*Mt*g0*cos(t1 - t3)*sin(t1 - t2)*sin(t3) + 5*m*r*t1*t1_d*cos(t1 - t2) + 2*m*r*t2_d*cos(t1 - t2)*sin(t1 - t3) - 4*Mt*l*t3*t3_d*sin(t1 - t2)*sin(t1 - t3) + 2*m*r*t2*t2_d*sin(t1 - t2)*sin(t1 - t3) - 4*Mt*r*t1*t1_d*cos(t1 - t2)*cos(t1 - t3)^2 - 4*Mt*r*t1*t1_d*cos(t1 - t3)*sin(t1 - t2)*sin(t1 - t3)))/(r*(4*Mh + 4*Mt + 5*m - 4*Mt*cos(t1 - t3)^2 - 4*m*cos(t1 - t2)^2)) - (2*(2*g0*m*cos(t1 - t2)*cos(t2) + 2*g0*m*sin(t1 - t2)*sin(t2) - m*r*t2_d*sin(t1 - t3) - 2*m*r*t1*t1_d*cos(t1 - t2)^2 + 2*m*r*t1*t1_d*sin(t1 - t2)^2))/(r*(4*Mh + 4*Mt + 5*m - 4*Mt*cos(t1 - t3)^2 - 4*m*cos(t1 - t2)^2)) - (16*m*cos(t1 - t2)*sin(t1 - t2)*(2*Mt*g0*cos(t1 - t3)*sin(t3) - 4*Mt*g0*sin(t1) - 3*g0*m*sin(t1) + 2*g0*m*cos(t1 - t2)*sin(t2) + 2*Mt*l*t3*t3_d*sin(t1 - t3) - m*r*t2*t2_d*sin(t1 - t3) + 2*Mt*r*t1*t1_d*cos(t1 - t3)*sin(t1 - t3) + 2*m*r*t1*t1_d*cos(t1 - t2)*sin(t1 - t2)))/(r*(- 4*m*cos(t1 - t2)^2 - 4*Mt*cos(t1 - t3)^2 + 4*Mh + 4*Mt + 5*m)^2) - (16*m*cos(t1 - t2)*sin(t1 - t2)*(5*g0*m*sin(t2) + 4*Mh*g0*sin(t2) + 4*Mt*g0*sin(t2) - 8*Mt*g0*cos(t1 - t2)*sin(t1) - 6*g0*m*cos(t1 - t2)*sin(t1) - 4*Mt*g0*cos(t1 - t3)^2*sin(t2) + 4*Mt*g0*cos(t1 - t2)*cos(t1 - t3)*sin(t3) + 4*Mh*r*t1*t1_d*sin(t1 - t2) + 4*Mt*r*t1*t1_d*sin(t1 - t2) + 5*m*r*t1*t1_d*sin(t1 - t2) + 4*Mt*l*t3*t3_d*cos(t1 - t2)*sin(t1 - t3) - 2*m*r*t2*t2_d*cos(t1 - t2)*sin(t1 - t3) - 4*Mt*r*t1*t1_d*cos(t1 - t3)^2*sin(t1 - t2) + 4*Mt*r*t1*t1_d*cos(t1 - t2)*cos(t1 - t3)*sin(t1 - t3)))/(r*(- 4*m*cos(t1 - t2)^2 - 4*Mt*cos(t1 - t3)^2 + 4*Mh + 4*Mt + 5*m)^2);
dLfh21 = (8*m*cos(t1 - t2)*sin(t1 - t2)*(5*g0*m*sin(t3) + 4*Mh*g0*sin(t3) + 4*Mt*g0*sin(t3) - 8*Mt*g0*cos(t1 - t3)*sin(t1) - 6*g0*m*cos(t1 - t3)*sin(t1) - 4*g0*m*cos(t1 - t2)^2*sin(t3) + 4*Mh*r*t1*t1_d*sin(t1 - t3) + 4*Mt*r*t1*t1_d*sin(t1 - t3) + 4*g0*m*cos(t1 - t2)*cos(t1 - t3)*sin(t2) + 5*m*r*t1*t1_d*sin(t1 - t3) + 4*Mt*l*t3*t3_d*cos(t1 - t3)*sin(t1 - t3) - 2*m*r*t2*t2_d*cos(t1 - t3)*sin(t1 - t3) - 4*m*r*t1*t1_d*cos(t1 - t2)^2*sin(t1 - t3) + 4*m*r*t1*t1_d*cos(t1 - t2)*cos(t1 - t3)*sin(t1 - t2)))/(l*(- 4*m*cos(t1 - t2)^2 - 4*Mt*cos(t1 - t3)^2 + 4*Mh + 4*Mt + 5*m)^2) - (8*g0*m*cos(t1 - t2)*sin(t1 - t2)*sin(t3) - 4*g0*m*cos(t1 - t2)*cos(t1 - t3)*cos(t2) - 4*g0*m*cos(t1 - t3)*sin(t1 - t2)*sin(t2) + 2*m*r*t2_d*cos(t1 - t3)*sin(t1 - t3) + 4*m*r*t1*t1_d*cos(t1 - t2)^2*cos(t1 - t3) - 4*m*r*t1*t1_d*cos(t1 - t3)*sin(t1 - t2)^2 + 8*m*r*t1*t1_d*cos(t1 - t2)*sin(t1 - t2)*sin(t1 - t3))/(l*(4*Mh + 4*Mt + 5*m - 4*Mt*cos(t1 - t3)^2 - 4*m*cos(t1 - t2)^2));
dLfh22 = -(2*(4*Mt*g0*cos(t1 - t2)*sin(t1 - t3)*sin(t3) - 8*Mt*g0*cos(t1 - t3)*sin(t1 - t3)*sin(t2) + 4*Mt*l*t3_d*cos(t1 - t2)*sin(t1 - t3) + 4*Mt*g0*cos(t1 - t2)*cos(t1 - t3)*cos(t3) - 4*Mt*l*t3*t3_d*cos(t1 - t2)*cos(t1 - t3) + 2*m*r*t2*t2_d*cos(t1 - t2)*cos(t1 - t3) - 4*Mt*r*t1*t1_d*cos(t1 - t2)*cos(t1 - t3)^2 + 4*Mt*r*t1*t1_d*cos(t1 - t2)*sin(t1 - t3)^2 - 8*Mt*r*t1*t1_d*cos(t1 - t3)*sin(t1 - t2)*sin(t1 - t3)))/(r*(4*Mh + 4*Mt + 5*m - 4*Mt*cos(t1 - t3)^2 - 4*m*cos(t1 - t2)^2)) - (2*(2*Mt*g0*sin(t1 - t3)*sin(t3) + 2*Mt*l*t3_d*sin(t1 - t3) + 2*Mt*g0*cos(t1 - t3)*cos(t3) - 2*Mt*l*t3*t3_d*cos(t1 - t3) + m*r*t2*t2_d*cos(t1 - t3) - 2*Mt*r*t1*t1_d*cos(t1 - t3)^2 + 2*Mt*r*t1*t1_d*sin(t1 - t3)^2))/(r*(4*Mh + 4*Mt + 5*m - 4*Mt*cos(t1 - t3)^2 - 4*m*cos(t1 - t2)^2)) - (16*Mt*cos(t1 - t3)*sin(t1 - t3)*(2*Mt*g0*cos(t1 - t3)*sin(t3) - 4*Mt*g0*sin(t1) - 3*g0*m*sin(t1) + 2*g0*m*cos(t1 - t2)*sin(t2) + 2*Mt*l*t3*t3_d*sin(t1 - t3) - m*r*t2*t2_d*sin(t1 - t3) + 2*Mt*r*t1*t1_d*cos(t1 - t3)*sin(t1 - t3) + 2*m*r*t1*t1_d*cos(t1 - t2)*sin(t1 - t2)))/(r*(- 4*m*cos(t1 - t2)^2 - 4*Mt*cos(t1 - t3)^2 + 4*Mh + 4*Mt + 5*m)^2) - (16*Mt*cos(t1 - t3)*sin(t1 - t3)*(5*g0*m*sin(t2) + 4*Mh*g0*sin(t2) + 4*Mt*g0*sin(t2) - 8*Mt*g0*cos(t1 - t2)*sin(t1) - 6*g0*m*cos(t1 - t2)*sin(t1) - 4*Mt*g0*cos(t1 - t3)^2*sin(t2) + 4*Mt*g0*cos(t1 - t2)*cos(t1 - t3)*sin(t3) + 4*Mh*r*t1*t1_d*sin(t1 - t2) + 4*Mt*r*t1*t1_d*sin(t1 - t2) + 5*m*r*t1*t1_d*sin(t1 - t2) + 4*Mt*l*t3*t3_d*cos(t1 - t2)*sin(t1 - t3) - 2*m*r*t2*t2_d*cos(t1 - t2)*sin(t1 - t3) - 4*Mt*r*t1*t1_d*cos(t1 - t3)^2*sin(t1 - t2) + 4*Mt*r*t1*t1_d*cos(t1 - t2)*cos(t1 - t3)*sin(t1 - t3)))/(r*(- 4*m*cos(t1 - t2)^2 - 4*Mt*cos(t1 - t3)^2 + 4*Mh + 4*Mt + 5*m)^2);
dLfh23 = (5*g0*m*cos(t3) + 4*Mh*g0*cos(t3) + 4*Mt*g0*cos(t3) - 8*Mt*g0*sin(t1 - t3)*sin(t1) - 6*g0*m*sin(t1 - t3)*sin(t1) - 4*g0*m*cos(t1 - t2)^2*cos(t3) - 4*Mh*r*t1*t1_d*cos(t1 - t3) - 4*Mt*r*t1*t1_d*cos(t1 - t3) - 5*m*r*t1*t1_d*cos(t1 - t3) + 4*g0*m*cos(t1 - t2)*sin(t1 - t3)*sin(t2) - 4*Mt*l*t3*t3_d*cos(t1 - t3)^2 + 4*Mt*l*t3*t3_d*sin(t1 - t3)^2 + 2*m*r*t2*t2_d*cos(t1 - t3)^2 - 2*m*r*t2*t2_d*sin(t1 - t3)^2 + 4*Mt*l*t3_d*cos(t1 - t3)*sin(t1 - t3) + 4*m*r*t1*t1_d*cos(t1 - t2)^2*cos(t1 - t3) + 4*m*r*t1*t1_d*cos(t1 - t2)*sin(t1 - t2)*sin(t1 - t3))/(l*(4*Mh + 4*Mt + 5*m - 4*Mt*cos(t1 - t3)^2 - 4*m*cos(t1 - t2)^2)) + (8*Mt*cos(t1 - t3)*sin(t1 - t3)*(5*g0*m*sin(t3) + 4*Mh*g0*sin(t3) + 4*Mt*g0*sin(t3) - 8*Mt*g0*cos(t1 - t3)*sin(t1) - 6*g0*m*cos(t1 - t3)*sin(t1) - 4*g0*m*cos(t1 - t2)^2*sin(t3) + 4*Mh*r*t1*t1_d*sin(t1 - t3) + 4*Mt*r*t1*t1_d*sin(t1 - t3) + 4*g0*m*cos(t1 - t2)*cos(t1 - t3)*sin(t2) + 5*m*r*t1*t1_d*sin(t1 - t3) + 4*Mt*l*t3*t3_d*cos(t1 - t3)*sin(t1 - t3) - 2*m*r*t2*t2_d*cos(t1 - t3)*sin(t1 - t3) - 4*m*r*t1*t1_d*cos(t1 - t2)^2*sin(t1 - t3) + 4*m*r*t1*t1_d*cos(t1 - t2)*cos(t1 - t3)*sin(t1 - t2)))/(l*(- 4*m*cos(t1 - t2)^2 - 4*Mt*cos(t1 - t3)^2 + 4*Mh + 4*Mt + 5*m)^2);

dLfh = [dLfh11, dLfh12, dLfh13;...
    dLfh21, dLfh22, dLfh23];  


K = [5, -1;...
    0, 12];

u = (dLfh*gx)\(-K*y - dLfh*fx);

dx(1:3) = x(4:6);
dx(4:6) = fx + gx*u;
dx = dx';

end

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