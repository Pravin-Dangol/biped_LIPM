function x_plus = impact_model_3link(x)

t_temp = x(1:3); t1 = t_temp(1); t2 = t_temp(2); t3 = t_temp(3);
td_temp = x(4:6); t1_d = td_temp(1); t2_d = td_temp(2); t3_d = td_temp(3);

[r,m,Mh,Mt,l,g0] = model_params_3link;

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