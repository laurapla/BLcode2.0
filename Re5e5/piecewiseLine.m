function y = piecewiseLine(x,a1,S1,S2)
% PIECEWISELINE   A line made of two pieces
% that is not continuous.

y = zeros(size(x));

% This example includes a for-loop and if statement
% purely for example purposes.
for i = 1:length(x)
    if abs(x(i)) < a1
        y(i) = 1-0.3*exp((abs(x(i))-a1)/S1);
    else
        y(i) = 0.04+0.66*exp(a1-(abs(x(i)))/S2);
    end
end
end