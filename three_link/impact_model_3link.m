function x_plus = impact_model_3link(x)

t_temp = num2cell(x(1:3)); [t1, t2, t3] = t_temp{:};
td_temp = num2cell(x(4:6)); [t1_d, t2_d, t3_d] = td_temp{:};

[r,m,Mh,Mt,l,~] = model_params_3link;

R = [0, 1, 0;...
    1, 0, 0;...
    0, 0, 1];

delq = R;

den = -3*m-4*Mh -2*Mt + 2*m*cos(2*t1 - 2*t2) + 2*Mt*cos(-2*t2 + 2*t3);

delq_d = zeros(3,3);
delq_d(1,1) = 1/den*(2*Mt*cos(-t1+2*t3-t2)-(2*m+4*Mh+2*Mt)*cos(t1-t2));
delq_d(1,2) = m/den;
delq_d(2,1) = 1/den*(m-(4*m+4*Mh+2*Mt)*cos(2*t1-2*t2)...
    +2*Mt*cos(2*t1-2*t3));
delq_d(2,2) = 1/den*2*m*cos(t1-t2);
delq_d(3,1) = r/(l*den)*(2*(m+Mh+Mt)*cos(t3+t1-2*t2)...
    -2*(m+Mh+Mt)*cos(-t1-t3)...
    +m*cos(-3*t1+2*t2+t3));
delq_d(3,2) = -r/(l*den)*cos(-t2+t3);
delq_d(3,3) = 1;

%Impact model:
x_plus = [delq*[t1; t2; t3];...
        delq_d*[t1_d; t2_d; t3_d]];

end