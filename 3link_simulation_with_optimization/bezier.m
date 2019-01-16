function b = bezier(s,M,a)

% Computes the Bezier polynomial

%for a normalized general coordinate s = (theta(q) - min(theta))/(max(theta) - min(theta))
%the function ouputs a 1 dimensional Bezier polynomial of order M
%the coefficient alpha of size 1 by M+1 must be provided

[i,j] = size(a);
if i > 1
    error('alpha vector must be 1D in length')
elseif j > M+1
    error('alpha vector must be size 1 by M+1')
end

%b = zeros(1,M+1); %b((1,k+1))
b = 0;
for k = 0:M
    b = b + a(1,k+1)*(factorial(M)/(factorial(k)*factorial(M-k)))*s^k*(1-s)^(M-k);
end

end