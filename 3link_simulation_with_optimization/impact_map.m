function [x_plus, f2] = impact_map(x)

%Provides impact map for a 3link, with q1 as cyclic variable

[r,m,Mh,Mt,l,~] = model_params_3link;

q = num2cell(x(1:3)); [q1, q2, q3] = q{:};

% De matrix
De=zeros(5,5);
De(1,1) = Mt*l^2 + Mh*r^2 + Mt*r^2 + (3*m*r^2)/2 - m*r^2*cos(q2) - 2*Mt*l*r*cos(q3);
De(1,2) = -(m*r^2*(2*cos(q2) - 1))/4;
De(1,3) = Mt*l*(l - r*cos(q3));
De(1,4) = Mt*l*cos(q1 + q3) - (3*m*r*cos(q1))/2 + (m*r*cos(q1 + q2))/2 - Mh*r*cos(q1) - Mt*r*cos(q1);
De(1,5) = Mt*l*sin(q1 + q3) - Mt*r*sin(q1) - (3*m*r*sin(q1))/2 - Mh*r*sin(q1) + (m*r*sin(q1 + q2))/2;
De(2,1) = De(1,2); 
De(2,2) = (m*r^2)/4; 
De(2,3) = 0;
De(2,4) = (m*r*cos(q1 + q2))/2;
De(2,5) = (m*r*sin(q1 + q2))/2;
De(3,1) = De(1,3); 
De(3,2) = De(2,3);
De(3,3) = Mt*l^2;
De(3,4) = Mt*l*cos(q1 + q3);
De(3,5) = Mt*l*sin(q1 + q3);
De(4,1) = De(1,4);
De(4,2) = De(2,4);
De(4,3) = De(3,4);
De(4,4) = Mh + Mt + 2*m;
De(5,1) = De(1,5);
De(5,2) = De(2,5);
De(5,3) = De(3,5);
De(5,5) = Mh + Mt + 2*m;

% E matrix
E=zeros(2,5);
E(1,1) = r*cos(q1 + q2) + r*cos(q1); 
E(1,2) = r*cos(q1 + q2); 
E(1,3) = 0;
E(1,4) = 1;
E(1,5) = 0;
E(2,1) = r*sin(q1 + q2) - r*sin(q1); 
E(2,2) = r*sin(q1 + q2); 
E(2,3) = 0;
E(2,4) = 0;
E(2,5) = 1;

%Impact map from pg56 (3.20)
Delta = [De -E';E zeros(2,2)]\[De*[x(4:6)';zeros(2,1)];zeros(2,1)]; %7x1

%Is this right? Check again
R = [1, 1, 0;...
    0, -1, 0;...
    0, -1, 1]; 
%q1 and q2 switch position after impact
x_plus(1:3) = (R*x(1:3)')';
x_plus(4:6) = (R*Delta(1:3))';

%x_plus(1) = x(2);x_plus(2) = x(1);x_plus(3) = x(3);
%x_plus(4) = temp(2);x_plus(5) = temp(1);x_plus(6) = temp(3);

f2 = Delta(5);
x_plus(7) = Delta(6);
x_plus(8) = Delta(7);
end