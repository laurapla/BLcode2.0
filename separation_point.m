function f = separation_point(alpha,alpha1,S1,S2)

% Function that computes the location of the separation point as a function
% of the angle of attack alpha
% f = separation point
% alpha = angle of attack in which we want to know the separation point
% [deg]
% alpha1 = alpha corresponding to f=0.7 [deg]
% S1 = coefficient that defines the stall characteristic [deg]
% S2 = coefficient that defines the stall characteristic [deg]

N = length(alpha);
f = zeros(1,N);

for i = 1:N
    if abs(alpha(i))<=alpha1
        f(i) = 1-0.3*exp((abs(alpha(i))-alpha1)/S1);
    else
        f(i) = 0.04+0.66*exp((alpha1-abs(alpha(i)))/S2);
    end
end

end