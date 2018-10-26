function [x_plus, f2] = impact_map_3link(x)

[r,m,Mh,Mt,l,~] = model_params_3link;

thetas = num2cell(x(1:3)); [t1, t2, t3] = thetas{:};

% De matrix
De=zeros(5,5);
De(1,1) = r^2*(5/4*m + Mh + Mt);
De(1,2) = -1/2*r^2*m*(cos(t1)*cos(t2) + sin(t1)*sin(t2));
De(1,3) = r*l*Mt*(cos(t3)*cos(t1) + sin(t3)*sin(t1));
De(1,4) = 1/2*r*cos(t1)*(3*m + 2*Mh + 2*Mt);
De(1,5) = -1/2*r*sin(t1)*(3*m + 2*Mh + 2*Mt);
De(2,1) = -1/2*r^2*m*(cos(t1)*cos(t2) + sin(t1)*sin(t2));
De(2,2) = 1/4*r^2*m;
De(2,4) = -1/2*m*r*cos(t2);
De(2,5) = 1/2*r*sin(t2)*m;
De(3,1) = r*l*Mt*(cos(t3)*cos(t1) + sin(t3)*sin(t1));
De(3,3) = l^2*Mt;
De(3,4) = l*cos(t3)*Mt;
De(3,5) = -Mt*l*sin(t3);
De(4,1) = 1/2*r*cos(t1)*(3*m + 2*Mh + 2*Mt);
De(4,2) = -1/2*m*r*cos(t2);
De(4,3) = l*cos(t3)*Mt;
De(4,4) = 2*m + Mh + Mt;
De(5,1) = -1/2*r*sin(t1)*(3*m + 2*Mh + 2*Mt);
De(5,2) = 1/2*r*sin(t2)*m;
De(5,3) = -Mt*l*sin(t3);
De(5,5) = 2*m + Mh + Mt;

% E matrix
E=zeros(2,5);
E(1,1) = r*cos(t1);
E(1,2) = -r*cos(t2);
E(1,4) = 1;
E(2,1) = -r*sin(t1);
E(2,2) = r*sin(t2);
E(2,5) = 1;

%Impact map from pg56 (3.21)
temp = [De -E';E zeros(2,2)]\[De*[x(4:6)';zeros(2,1)];zeros(2,1)]; %7x1

R = [0, 1, 0;...
    1, 0, 0;...
    0, 0, 1]; 
%theta1 and theta 2 switch position after impact
x_plus(1) = x(2);
x_plus(2) = x(1);
x_plus(3) = x(3);
x_plus(4) = temp(2);
x_plus(5) = temp(1);
x_plus(6) = temp(3);
x_plus(7) = temp(6);
x_plus(8) = temp(7);
f2 = temp(5);

end