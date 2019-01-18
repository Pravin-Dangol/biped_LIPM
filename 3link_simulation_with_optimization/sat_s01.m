function y = sat_s01(x)

if x > 1
    y = 1;
elseif x < 0
    y = 0;
else
    y = x;
end

end