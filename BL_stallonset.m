function Cnprime = BL_stallonset(s,T_P,Cnc,Cnp)

% Function that computes the equivalent critical normal force coefficient
% for an unsteady motion
% Cnprime = equivalent critical normal force coefficient
% s = non-dimensional time vector
% T_P = time constant for leading edge pressure response
% Cnc = normal circulatory force coefficient (attached flow)
% Cnp = normal force coefficient (attached flow)

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