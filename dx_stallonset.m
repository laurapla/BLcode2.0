function dx = dx_stallonset(s,x,T_P,Cnp,s_time)

% Function that gives the differential equation of the critical pressure
% (normal force coefficient) criterion
% s = non-dimensional time
% T_P = time constant for leading edge pressure response [s]
% Cnp = vector of normal force coefficients of the attached flow
% s_time = non-dimensional time steps at which Cnp was calculated

CNp = interp1(s_time,Cnp,s);

dx = -x/T_P+CNp/T_P;

end