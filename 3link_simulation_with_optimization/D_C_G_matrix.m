function [D,C,G,B] = state_matrix(x)

[r,m,Mh,Mt,l,g] = model_params_3link;

q_temp = num2cell(x(1:3)); [q1, q2, q3] = q_temp{:};
dq_temp = num2cell(x(4:6)); [dq1, dq2, dq3] = dq_temp{:};

%Mass matrix
D = zeros(3,3);
D(1,1) = Mt*l^2 + Mh*r^2 + Mt*r^2 + (3*m*r^2)/2 - m*r^2*cos(q2) - 2*Mt*l*r*cos(q3);
D(1,2) = -(m*r^2*(2*cos(q2) - 1))/4;
D(1,3) = Mt*l*(l - r*cos(q3));
D(2,1) = D(1,2); 
D(2,2) = (m*r^2)/4; 
D(2,3) = 0;
D(3,1) = D(1,3); 
D(3,2) = D(2,3);
D(3,3) = Mt*l^2; 

%C matrix
C = zeros(3,3);
C(1,1) = (dq2*m*sin(q2)*r^2)/2 + Mt*dq3*l*sin(q3)*r;
C(1,2) = (m*r^2*sin(q2)*(dq1 + dq2))/2;
C(1,3) = Mt*l*r*sin(q3)*(dq1 + dq3);
C(2,1) = -(dq1*m*r^2*sin(q2))/2;
C(2,2) = 0; 
C(2,3) = 0;
C(3,1) = -Mt*dq1*l*r*sin(q3);
C(3,2) = 0; 
C(3,3) = 0;

%G matrix
G = [-(g*(2*Mh*r*sin(q1) + 2*Mt*r*sin(q1) + 3*m*r*sin(q1) - 2*Mt*l*sin(q1 + q3) - m*r*sin(q1 + q2)))/2;...
    (g*m*r*sin(q1 + q2))/2;...
    Mt*g*l*sin(q1 + q3)];

%B matrix
B = [0, 0;...
    1, 0;...
    0, 1];

end