function fm = BL_mseparation(s_span,alpha,T_f,alpha1,S1,S2)

% Time vector (to solve the dynamics)
s_time = s_span;

% Initial condition
xi = 1;

% System of equations
[~,x] = ode45(@(t,x) dx_separation(t,x,T_f,alpha1,S1,S2,alpha,s_time),s_span,xi);

% Output
fm = x.';

end