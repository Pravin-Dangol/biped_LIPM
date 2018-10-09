function x_plus = impack_model_3link(x)

[q1, q2, q3] = x(1:3);

%setting parameters from biped book
m = 5; Mh = 15; Mt = 10;
r = 1; l = 0.5; g0 = 9.81;

%conversion between theta and q

q2t = [1, 0, 1;...
    0, 1, 1;...
    0, 0, 1];

T1_3 = q2t*[q1; q2; q3] + [pi; pi; 0];

t1 = T1_3(1);
t2 = T1_3(2);
t3 = T1_3(3);

R = [0, 1, 0;...
    1, 0, 0;...
    0, 0, 1];

den = -3*m-4*Mh -2*Mt + 2*m*cos(2*t1 - 2*t2) + 2*Mt*cos(-2*t2 + 2*t3);

dqd11 = 1/den*(2*Mt*cos(-t1+2*t3-t2)-(2*m+4*Mh+2*Mt)*cos(t1-t2));
dqd12 = m/den;
dqd13 = 0;
dqd21 = 1/den*(m-(4*m+4*Mh+2*Mt)*cos(2*t1-2*t2)...
    +2*Mt*cos(2*t1-2*t3));
dqd22 = 1/den*2*m*cos(t1-t2);
dqd23 = 0;
dqd31 = r/(l*den)*(2*(m+Mh+Mt)*cos(t3+t1-2*t2)...
    -2*(m+Mh+Mt)*cos(-t1-t3)...
    +m*cos(-3*t1+2*t2+t3));
dqd32 = -r/(l*den)*cos(-t2+t3);
dqd33 = 1;

delq = R;

delq_d = [dqd11, dqd12, dqd13;...
    dqd21, dqd22, dqd23;...
    dqd31, dqd32, dqd33];

%Impact model:

x_plus = [delq*[t1; t2; t3];...
        delq_d*[t1_d; t2_d; t3_d]];

end