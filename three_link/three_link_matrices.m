function [D,C,G,B,fx,gx] = three_link_matrices(x)

[r,m,Mh,Mt,l,g0] = model_params_3link;

t_temp = x(1:3); t1 = t_temp(1); t2 = t_temp(2); t3 = t_temp(3);
td_temp = x(4:6); t1_d = td_temp(1); t2_d = td_temp(2); t3_d = td_temp(3);

%Mass matrix

D11 = (5/4*m + Mh + Mt)*r^2;
D12 = -0.5*m*r^2*cos(t1-t2);
D13 = Mt*r*l*cos(t1-t3);
D22 = 0.25*m*r^2; D21 = D12; D23 = 0;
D33 = Mt*l^2; D31 = D13; D32 = D23;

D = [D11, D12, D13;...
    D21, D22, D23;...
    D31, D32, D33];

%C matrix

C11 = 0;
C12 = -0.5*m*r^2*sin(t1-t3)*t2_d;
C13 = Mt*r*l*sin(t1-t3)*t3_d;
C21 = 0.5*m*r^2*sin(t1-t2)*t1_d;
C22 = 0; C23 = 0;
C31 = -Mt*r*l*sin(t1-t3)*t1_d;
C32 = 0; C33 = 0;

C = [C11, C12, C13;...
    C21, C22, C23;...
    C31, C32, C33];

%G matrix

G = [-0.5*g0*(2*Mt+3*m+2*Mt)*r*sin(t1);...
    0.5*g0*m*r*sin(t2);...
    -g0*Mt*l*sin(t3)];

%B matrix

B = [-1, 0;...
    0, -1;...
    1, 1];

fx = D\(-C*x(4:6)-G);
gx = D\B;

end