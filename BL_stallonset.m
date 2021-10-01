function Cnprime = BL_stallonset(s,T_P,Cnc,Cnp)

% Function that computes the equivalent critical normal force coefficient
% for an unsteady motion
% Cnprime = equivalent critical normal force coefficient
% s_span = non-dimensional time vector
% T_P = time constant for leading edge pressure response
% C_N1 = critical normal force coefficient
% Cnp = normal force coefficient for attached flow

% % Initial condition
% xi = 0;
% 
% % Non-dimensional time vector (to solve the dynamics)
% s_time = s;
% 
% % System of equations
% [~,x] = ode45(@(t,x) dx_stallonset(t,x,T_P,Cnp,s_time),s,xi);
% 
% % Output
% Cnprime = x.';

N = size(s,2);

Dp = zeros(1,N);
Cnprime = zeros(1,N);

for i = 2:N
    
    ds = s(i)-s(i-1);
    Ep = exp(-ds/T_P);
    Dp(i) = Dp(i-1)*Ep+(Cnc(i)-Cnc(i-1))*Ep^0.5;
    Cnprime(i) = Cnp(i)-Dp(i);
    
end


end