function dx = dx_TEseparation(t,x,fprime,T_f,s_time)

% Function that gives the differential equation to compute the effects of
% the trailing edge separation
% fprime = effective separation point
% T_f = time constant for separation point movement [s]
% s_time = non-dimensional time instants at which fprime was computed

ffprime = interp1(s_time,fprime,t);

dx = -x/T_f+ffprime/T_f;

end