function alpha_E = effective_angle(beta,V,c,x1,x2)

% Function that computes the effective angle of attack of the airfoil due
% to the shed wake (circulatory) terms
% alpha_E = effective angle of attack of the airfoil due to the shed wake
% beta = compressibility factor
% V = free-stream velocity
% c = airfoil chord
% x1 = first state of the system
% x2 = second state of the system

A1 = 0.3; A2 = 0.7;
b1 = 0.14; b2 = 0.53;

alpha_E = beta^2*(2*V/c)*(A1*b1*x1+A2*b2*x2);

end