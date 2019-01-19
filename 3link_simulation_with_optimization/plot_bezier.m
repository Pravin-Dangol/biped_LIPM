
function plot_bezier(a)

% Plot Bezier polynomial based on given coefficient

s = 0:0.01:1;

b = a(1)*(s-1).^4 - 4*a(2)*s.*(s-1).^3 + 6*a(3)*s.^2.*(s-1).^2 - 4*a(4)*s.^3.*(s-1) + a(5)*s.^4 ;

plot(s,b)

hold on

plot(0:1/4:1,a,'o')

end