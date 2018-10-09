function [D,C,G,B,fx,gx] = three_link_matrices(x(1:6))

%setting parameters from biped book
m = 5; Mh = 15; Mt = 10;
r = 1; l = 0.5; g0 = 9.81;

%Mass matrix

D11 = (5/4*m + Mh + Mt)*r^2;
D12 = -0.5*m*r^2*cos(q1-q2);
D13 = (5/4*m + Mh + Mt -m/2*cos(q1-q2))*r^2 - Mt*r*l*cos(q1);
D21 = D12;
D22 = 0.25*m*r^2; 
D23 = (m/4 - m/2*cos(q1-q2))*r^2;
D31 = D13; 
D32 = D23;
D33 = (Mh + 3/2*m + Mt -m*cos(q1-q2))*r^2 - 2*Mt*r*l*cos(q1) + Mt*l^2; 


D = [D11, D12, D13;...
    D21, D22, D23;...
    D31, D32, D33];

%C matrix

C11 = 0;
C12 = -0.5*m*r^2*sin(q1-q2)*(q2_d+q3_d);
C13 = -r/2*(m*r*sin(q1-q2)*q2_d + m*r*sin(q1-q2)*q3_d + 2*Mt*l*sin(q1)*q3_d);
C21 = 0.5*m*r^2*sin(q1-q2)*(q1_d+q3_d);
C22 = 0; 
C23 = C21;
C31 = 0.5*(m*r^2*sin(q1-q2)+ 2*Mt*r*l*sin(q1))*(q1_d+q3_d);
C32 = -0.5*m*r^2*sin(q1-q2)*(q2_d+q3_d); 
C33 = 0.5*(m*r^2*sin(q1-q2)*q1_d + 2*Mt*r*l*sin(q1)*q1_d - m*r^2*sin(q1-q2)*q2_d);

C = [C11, C12, C13;...
    C21, C22, C23;...
    C31, C32, C33];

%G matrix

G1 = 0.5*g0*(3*m + 2*Mh + 2*Mt)*r*sin(q1);
G2 = -0.5*g0*m*r*sin(q2+q3);
G3 = 0.5*g0*((3*m + 2*Mh + 2*Mt)*r*sin(q1+q3) - m*r*sin(q1+q3)) - g0*Mt*l*sin(q3);

G = [G1;...
    G2;...
    G3];

%B matrix

B = [-1, 0;...
    0, -1;...
    0, 0]; 

fx = inv(D)*(-C*x(4:6)-G);
gx = inv(D)*B;

end