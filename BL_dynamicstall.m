function [Cnv, Cmv] = BL_dynamicstall(s_span,fprimeprime,alpha,alpha_E,Cnprime,C_Nalpha,T_v,T_vl,C_N1)

% Function that computes the contribution of the vortex to the normal force
% coefficient
% Cnv = normal force coefficient due to leading-edge vortex (vector)
% Cmv = moment coefficient due to leading-edge vortex (vector)
% s_span = non-dimensional time (vector)
% fprimeprime = unsteady separation point (vector)
% alpha = angle of attack (vector) [rad]
% alpha_E = effective angle of attack (vector) [rad]
% Cnprime = equivalent critical normal force coefficient (vector)
% C_Nalpha = normal force surve slope [1/rad]
% T_v = time constant for vortex lift
% T_vl = time constant for vortex lift traverse
% C_N1 = critical normal force coefficient

% Time rate of change of circulation
tauv = vortex_time(Cnprime,s_span,C_N1);

Cnc = C_Nalpha*alpha_E; % Circulatory normal force coefficient
Kn = (1+sqrt(fprimeprime)).^2/4;

N = size(s_span,2);

Cv = zeros(1,N);
Cnv = zeros(1,N);
Ds = zeros(1,N);
Sa = zeros(1,N);
df = zeros(1,N);

for i = 2:N
    
    ds_span = s_span(i)-s_span(i-1);
    Sa(i) = alpha(i)-alpha(i-1);
    df(i) = fprimeprime(i)-fprimeprime(i-1);
    sigma = sigmav(tauv(i), T_vl, Sa(i), df(i));
    
    if abs(Cnprime(i))<C_N1 && df(i)<0
        Ds(i) = 1;
    elseif tauv(i)>0 && tauv(i)<=T_vl
        Ds(i) = 1;
    elseif Sa(i)>0 && df(i)>0
        Ds(i) = 1;
    end
    
    Cv(i) = Ds(i)*Cnc(i)*(1-Kn(i));
    Ev = exp(-sigma*ds_span/T_v);
    Cnv(i) = Cnv(i-1)*Ev+Ds(i)*(Cv(i)-Cv(i-1))*Ev^0.5;
    
end

n_t = 100;
Ds2 = ones(1,100);
for i = N-n_t+1:N
    if Cv(i-1)==0 && Cv(i-2)==0
        Ds2(i-N+n_t) = 0;
    elseif i<N-2
        if Cv(i+1)==0 && Cv(i+2)==0
            Ds2(i-N+n_t) = 0;
        end
    end
end
% dCv = Ds2.*gradient(Cv(N-n_t+1:N))./gradient(s_span(N-n_t+1:N));

% Center of pressure
CPv = 0.25*(1-cos(pi*tauv/T_vl));

% Pitching moment coefficient
Cmv = -CPv.*Cnv;

end