function dx = dx_separation(t,x,T_f,alpha1,S1,S2,alpha,time)

% Function that gives the differential equation that computes the effective
% separation point to compute the mean center of pressure during flow
% reattachment
% T_f = time constant for separation point movement [s]
% alpha1 = angle of attack at which the separation poInt is at f=0.7 [rad]
% S1 = stall characteristic coefficient for alpha<alpha1
% S2 = stall characteristic coefficient for alpha>alpha1
% alpha = angle of attack [rad]
% time = time vector at which alpha_ef was computed [s]

% Separation point
f = separation_point(alpha,alpha1,S1,S2);
fQS = interp1(time,f,t);

dx = -2*x/T_f+2*fQS/T_f;

end