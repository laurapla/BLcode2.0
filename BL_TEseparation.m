function [Cnf, Cmf, Cc, fprimeprime, fprime] = BL_TEseparation(s_span,Cnprime,C_Nalpha,C_N1,alpha0,alpha,alpha_E,alpha1,dalpha1,S1,S2,T_f,T_vl,K0,K1,K2,eta,D_f)

% Function that computes the variations in the normal force, pitching
% moment and drag force coefficients due to trailing edge separation
% Cnf = normal force coefficient due to trailing edge separation
% Cmf = pitching moment due to trailing edge separation
% Cdf = drag force due to trailing edge separation
% fprimeprime = unsteady trailing edge separation point
% s_span = non-dimensional time vector
% Cnprime = equivalent critical normal force coefficient (vector)
% C_Nalpha = normal force coefficient curve slope
% Cm0 = zero lift pitching moment coefficient
% alpha = angle of attack (vector) [rad]
% alpha_E = effective angle of attack (vector) [rad]
% alpha1 = angle of attack at which the separation point is f=0.7
% S1 = coefficient that defines the stall characteristic
% S2 = coefficient that defines the stall characteristic
% T_f = time constant for separation point movement [s]
% K0 = aerodynamic center offset from the 1/4-chord
% K1 = effect on the center of pressure due to the growth of separated flow
% K2 = describes the shape of the moment break at stall
% eta = viscous effects factor

% Effective angle of attack (for separation)
alpha_f = Cnprime./C_Nalpha+alpha0;

tauv = vortex_time(Cnprime,s_span,C_N1);

tolerance = 2e-4;

N = size(s_span,2);
Df = zeros(1,N);
fprime = ones(1,N);
fprimeprime = ones(1,N);
Dfr = zeros(1,N);
fM = ones(1,N);
fr = ones(1,N);

for i = 2:N
    
    % Trailing edge separation point
    ds = s_span(i)-s_span(i-1);
    Sa = alpha(i)-alpha(i-1);
    if Sa<0
        Da1 = (1-fprimeprime(i-1))^0.25*dalpha1;
    else
        Da1 = 0;
    end
    aalpha1 = alpha1-Da1;
    fprime(i) = separation_point(alpha_f(i),aalpha1,S1,S2);
    
    % Pre-calculation for reattachment separation point
    if Sa>=0
        fM(i) = fprime(i);
    else
        fM(i) = separation_point(alpha(i),aalpha1,S1,S2);
    end
    
    fprimeprime_ant = 0.8*fprimeprime(i-1); fr_ant = 0.8*fr(i-1); resta = 100;
    while resta>tolerance
        
        [sigma1, sigma3] = sigma13(Cnprime(i), C_N1, Sa, fprimeprime_ant-fprimeprime(i-1), fprimeprime_ant, fr_ant, tauv(i), T_vl);
        
        % Unsteady trailing edge separation point (unsteady boundary layer
        % effects)
        Ef = exp(-sigma1*ds/T_f);
        Df(i) = Df(i-1)*Ef+(fprime(i)-fprime(i-1))*Ef^0.5;
        fprimeprime(i) = fprime(i)-Df(i);
        
        % Reattachment separation point
        Em = exp(-sigma3*ds/T_f);
        Dfr(i) = Dfr(i-1)*Em+(fM(i)-fM(i-1))*Em^0.5;
        fr(i) = fM(i)-Dfr(i);
        
        resta = abs(fprimeprime(i)-fprimeprime_ant);
        fprimeprime_ant = fprimeprime(i);
        fr_ant = fr(i);
        
    end
    
end

% Normal force coefficient
Cnc = C_Nalpha*alpha_E; % Circulatory normal force coefficient
Cnf = Cnc.*((1+sqrt(fprimeprime))/2).^2;

% Pitching moment coefficient
m = 2;
Cmf = (K0+K1*(1-fr)+K2*sin(pi*fr.^m)).*Cnc.*(1+sqrt(fr)).^2/4;

% Drag force coefficient
Phi = zeros(1,N);
for i = 1:N
    if abs(Cnprime(i))<=C_N1
        Phi(i) = 1;
    else
        exponent = min(D_f*(abs(Cnprime(i))-C_N1),1);
        Phi(i) = fprimeprime(i)^exponent;
    end
end
Cc = eta*Cnc.*tan(alpha_E+alpha0).*sqrt(fprimeprime).*Phi; % Chord force coefficient

end