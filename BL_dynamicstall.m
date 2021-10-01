function [Cnv, Cmv, Cdv] = BL_dynamicstall(s_span,fprimeprime,alpha,alpha_E,Cnprime,C_Nalpha,T_v,T_vl,C_N1)

% Function that computes the contribution of the vortex that is 

% Time rate of change of circulation
tauv = vortex_time(Cnprime,s_span,C_N1);

Cnc = C_Nalpha*alpha_E; % Circulatory normal force coefficient
Kn = (1+sqrt(fprimeprime)).^2/4;

N = size(s_span,2);

Cv = zeros(1,N);
Cnv = zeros(1,N);

for i = 2:N
    
    ds_span = s_span(i)-s_span(i-1);
    Sa = alpha(i)-alpha(i-1);
    df = fprimeprime(i)-fprimeprime(i-1);
    sigma = sigma2(tauv(i), T_vl, Sa, df);
    
    if tauv(i)>0 && tauv(i)<=T_vl
        Ds = 1;
    elseif abs(Cnprime(i))<C_N1 && df<0
        Ds = 1;
    elseif Sa>0 && df>0
        Ds = 1;
    else
        Ds = 0;
    end
    
    Cv(i) = Ds*Cnc(i)*(1-Kn(i));
    Ev = exp(-sigma*ds_span/T_v);
    Cnv(i) = Cnv(i-1)*Ev+Ds*(Cv(i)-Cv(i-1))*Ev^0.5;
    
end

% Center of pressure
CPv = 0.25*(1-cos(pi*tauv/T_vl));

% Pitching moment coefficient
Cmv = -CPv.*Cnv;

% Drag force coefficient
Cdv = ones(1,N);

end