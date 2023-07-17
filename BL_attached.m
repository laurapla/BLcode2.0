function [Cnp, Cmp, Ccp, alpha_E, Cni] = BL_attached(t_span,alpha,q,V,M,c,e,C_Nalpha,alpha0,Cm0,x_ac)

% Function that computes the normal force coefficient, the moment
% coefficient and the drag force coefficient of the attached flow
% Cnp = normal force coefficient (vector)
% Cmp = moment coefficient (vector)
% Ccp = chord force coefficient (vector)
% alpha_E = effective angle of attack (Cnc/C_Nalpha) [rad]
% Cn_i = normal non-circulatory force coefficient (vector)
% t_span = times at which the outputs will be computed (vector) [s]
% alpha = angle of attack (vector) [rad]
% q = pitching rate (vector) [rad/s]
% V = free-stream velocity [m/s]
% M = Mach number
% c = airfoil chord [m]
% e = pitching axis normalized by the chord
% C_Nalpha = normal force curve slope [1/rad]
% alpha0 = angle of attack of zero lift [rad]
% Cm0 = zero lift pitching moment coefficient
% x_ac = aerodynamic center normalized by the chord

% Prandtl-Glauert compressibility correction factor
beta = sqrt(1-M^2);

% Constants
A1 = 0.636; A2 = 1-A1; A3 = 1.5; A4 = 1-A3; A5 = 1.0;
b1 = 0.339; b2 = 0.249; b3 = 0.25; b4 = 0.1; b5 = 0.5;
kNa = 0.75; kNq = 0.75; kMa = 0.8; kMq = 0.8;

a = V/M; % Speed of sound
T_I = c/a; % Basic non-circulatory time constant

% Constants
K_alpha = kNa/((1-M)+0.5*C_Nalpha*beta^2*M^2*(A1*b1+A2*b2));
K_q = kNq/((1-M)+C_Nalpha*beta^2*M^2*(A1*b1+A2*b2));
K_alphaM = kMa*(A3*b4+A4*b3)/(b3*b4*(1-M));
K_qM = 7*kMq/(15*(1-M)+1.5*C_Nalpha*beta^2*M^2*b5);

T_alpha = K_alpha*T_I;
T_q = K_q*T_I;

% Dimensions
s_span = 2*V*t_span/c;
N = size(s_span,2);
dat = zeros(1,N);
dqt = zeros(1,N);
X1 = zeros(1,N);
X2 = zeros(1,N);
X3 = zeros(1,N);
X4 = zeros(1,N);
X5 = zeros(1,N);
X6 = zeros(1,N);
X7 = zeros(1,N);
X8 = zeros(1,N);



% NORMAL FORCE COEFFICIENT
for i = 2:N
    
    ds = s_span(i)-s_span(i-1);
    dt = t_span(i)-t_span(i-1);
    da = alpha(i)-alpha(i-1);
    dq = q(i)-q(i-1);
    dat(i) = da/dt;
    dqt(i) = dq/dt;
    
    % Circulatory - AOA + pitch rate contribution
    % The input of the ODEs of states 1 and 2 is the angle of attack at the
    % 3/4-chord (Leishman, J.G. Principles of Helicopter Aerodynamics)
    X1(i) = X1(i-1)*exp(-b1*beta^2*ds)+A1*exp(-b1*beta^2*ds/2)*(da+dq*(1/4-e));
    X2(i) = X2(i-1)*exp(-b2*beta^2*ds)+A2*exp(-b2*beta^2*ds/2)*(da+dq*(1/4-e));
    
    % Non-circulatory - AOA contribution
    X3(i) = X3(i-1)*exp(-dt/T_alpha)+(dat(i)-dat(i-1))*exp(-dt/(2*T_alpha));
    
    % Non-circulatory - Pitch rate contribution
    X4(i) = X4(i-1)*exp(-dt/T_q)+(dqt(i)-dqt(i-1))*exp(-dt/(2*T_q));
    
end

% Circulatory normal force
alpha_E= alpha+q*(1/4-e)-alpha0-X1-X2;
Cnc = C_Nalpha*alpha_E;

% Non-circulatory normal force
Cni_a = 4*T_alpha/M*(dat-X3);
Cni_q = -T_q/M*(dqt-X4);
Cni = Cni_a+Cni_q;

% Total normal force under attached conditions
Cnp = Cnc+Cni;


% AERODYNAMIC MOMENT COEFFICIENT
for i = 2:N
    
    ds = s_span(i)-s_span(i-1);
%     da = alpha(i)-alpha(i-1);
    dq = q(i)-q(i-1);
    dt = t_span(i)-t_span(i-1);

    % Non-circulatory - AOA contribution
    X5(i) = X5(i-1)*exp(-dt/(b3*K_alphaM*T_I))+A3*da*exp(-dt/(2*b3*K_alphaM^2*T_I));
    X6(i) = X6(i-1)*exp(-dt/(b4*K_alphaM*T_I))+A4*da*exp(-dt/(2*b4*K_alphaM^2*T_I));

    % Circulatory - Pitch rate contribution
    X7(i) = X7(i-1)*exp(-b5*beta^2*ds)+A5*dq*exp(-b5*beta^2*ds/2);
    
    % Non-circulatory - Pitch rate contribution
    X8(i) = X8(i-1)*exp(-dt/(K_qM^2*T_I))+(dqt(i)-dqt(i-1))*exp(-dt/(2*K_qM^2*T_I));
    
end

% Circulatory moment
Cmc_a = -Cnc*(x_ac-0.25);
Cmc_q = -C_Nalpha*b5*beta^2/16*X7;
Cmc = Cmc_a+Cmc_q;

% Non-circulatory moment
Cmi_a = -Cni/4;
Cmi_q = -7*K_qM^2*T_I/(12*M)*(dqt-X8);
Cmi = Cmi_a+Cmi_q;

% Total aerodynamic moment coefficient
Cmp = Cmc+Cmi+Cm0;


% CHORD FORCE
Ccp = Cnc.*tan(alpha_E+alpha0);

end