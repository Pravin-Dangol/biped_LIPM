function [D,C,G,B] = state_matrix_3link(t,x)

[r,m,Mh,Mt,l,g0] = model_params_3link;

t_temp = num2cell(x(1:3)); [t1, t2, t3] = t_temp{:};
td_temp = num2cell(x(4:6)); [t1_d, t2_d, t3_d] = td_temp{:};

%Mass matrix
D = zeros(3,3);
D(1,1) = (5/4*m + Mh + Mt)*r^2;
D(1,2) = -1/2*m*r^2*cos(t1-t2);
D(1,3) = Mt*r*l*cos(t1-t3);
D(2,1) = D(1,2); 
D(2,2) = 1/4*m*r^2; 
D(2,3) = 0;
D(3,1) = D(1,3); 
D(3,2) = D(2,3);
D(3,3) = Mt*l^2; 

%C matrix
C = zeros(3,3);
C(1,1) = 0;
C(1,2) = -1/2*m*r^2*sin(t1-t3)*t2_d;
C(1,3) = Mt*r*l*sin(t1-t3)*t3_d;
C(2,1) = 1/2*m*r^2*sin(t1-t2)*t1_d;
C(2,2) = 0; 
C(2,3) = 0;
C(3,1) = -Mt*r*l*sin(t1-t3)*t1_d;
C(3,2) = 0; 
C(3,3) = 0;

%G matrix
G = [-1/2*g0*(2*Mt+3*m+2*Mt)*r*sin(t1);...
    1/2*g0*m*r*sin(t2);...
    -g0*Mt*l*sin(t3)];
    
%B matrix
B = [-1, 0;...
    0, -1;...
    1, 1];

end
    