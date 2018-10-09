% bezier polynomial of degree 2
%
% by Alireza Ramezani, 10-6-2018, Pawtucket, RI
function [y,p1,p2,p3] = func_bezier(stp_ln,stp_ht,n)


%
s = linspace(0,stp_ln,n);

% equally placed points 
p1 = [0,0].';
p2 = [1/2,stp_ht].';
p3 = [1,0].';

y = [];

for i=1:n
    tmp = (1-s(i)/stp_ln)^2*p1+2*(1-s(i)/stp_ln)*s(i)/stp_ln*p2+s(i)^2/stp_ln^2*p3;
    y = [y;[s(i),tmp(2)]];
end

end
