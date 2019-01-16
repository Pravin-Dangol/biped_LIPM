function q = map_z_to_x(z,a)

q1 = z(1); dq1 = z(2); 
a2 = a(1:5); a3 = a(6:end);

a21 = a(1); a22 = a(2); a23 = a(3); a24 = a(4); a25 = a(5); 
a31 = a(6); a32 = a(7); a33 = a(8); a34 = a(9); a35 = a(10); 

M = 4;

delta_theta = deg2rad(30);              %difference between min and max q1
s = (q1 + delta_theta/2)/delta_theta;   %normalized general coordinate

q2 = bezier(s,M,a2);
q3 = bezier(s,M,a3)';

dq2 = d_ds_bezier(s,M,a2)*dq1/delta_theta;
dq3 = d_ds_bezier(s,M,a3)*dq1/delta_theta;

q = [z(1), q2, q3, z(2), dq2, dq3];

end