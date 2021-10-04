function [Cnp, Cmp, Ccp, alpha_E, Cni] = BL_attached(t_span,alpha,q,V,M,c,C_Nalpha,alpha0,x_ac)

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
% C_Nalpha = normal force curve slope [1/rad]
% % x_ac = aerodynamic center normalized by the chord
% 
% % Time vector (to solve the dynamics)
% time = t_span;
% 
% % Matrices of the system
% [A,B,C,D] = attached_matrices(V,M,c,C_Nalpha,x_ac);
% 
% % Inputs of the system
% u = [alpha; q];
% 
% % Initial conditions
% xi = zeros(8,1);
% 
% % System of equations
% opts = odeset('RelTol',1e-5,'AbsTol',1e-8);
% [~,x] = ode45(@(t,x) dx_attached(t,x,A,B,alpha,q,time),t_span,xi,opts);
% x = x.';
% 
% % Output
% y = C*x+D*u;
% Cnp = y(1,:);
% Cmp = y(2,:);
% 
% % Circulatory and non-circulatory components
% alpha_E = effective_angle(sqrt(1-M^2),V,c,x(1,:),x(2,:)); % Effective angle [rad]
% Cn_i = Cnp-C_Nalpha*(alpha_E);
% 
% % Drag force
% Ccp = C_Nalpha*sin(alpha_E).^2; % Chord force coefficient

% Prandtl-Glauert compressibility correction factor
beta = sqrt(1-M^2);

% Constants
A1 = 0.636; A2 = 1-A1; A3 = 1.5; A4 = -0.5; A5 = 1.0;
b1 = 0.339; b2 = 0.249; b3 = 0.25; b4 = 0.1; b5 = 0.5;
kNa = 0.75; kNq = 0.75; kMa = 0.8; kMq = 0.8;

a = V/M; % Speed of sound
T_I = c/a; % Basic non-circulatory time constant

% Constants
K_alpha = kNa/((1-M)+C_Nalpha*beta*M^2*(A1*b1+A2*b2));
K_q = kNq/(0.5*(1-M)+2*pi*beta*M^2*(A1*b1+A2*b2));
K_qM = 7*kMq/(15*(1-M)+3*pi*beta*M^2*A5*b5);

T_alpha = K_alpha*T_I;
T_q = K_q*T_I;

C_Nalphac = C_Nalpha/beta;

% Dimensions
s_span = 2*V*t_span/c;
N = size(s_span,2);
X1 = zeros(1,N);
X2 = zeros(1,N);
X3 = zeros(1,N);
X4 = zeros(1,N);
dat = zeros(1,N);
dqt = zeros(1,N);
Kprime_a = zeros(1,N);
Kprime_q = zeros(1,N);
Kprime_qM = zeros(1,N);
Kprimeprime_qM = zeros(1,N);


% NORMAL FORCE COEFFICIENT
for i = 2:N
    
    ds = s_span(i)-s_span(i-1);
    dt = t_span(i)-t_span(i-1);
    da = alpha(i)-alpha(i-1);
    dq = q(i)-q(i-1);
    dat(i) = da/dt;
    dqt(i) = dq/dt;
    
    % Circulatory - AOA contribution
    X1(i) = X1(i-1)*exp(-b1*beta^2*ds)+A1*exp(-b1*beta^2*ds/2)*da;
    X2(i) = X2(i-1)*exp(-b2*beta^2*ds)+A2*exp(-b2*beta^2*ds/2)*da;
    
    % Circulatory - Pitch rate contribution
    X3(i) = X3(i-1)*exp(-b1*beta^2*ds)+A1*exp(-b1*beta^2*ds/2)*dq;
    X4(i) = X4(i-1)*exp(-b2*beta^2*ds)+A2*exp(-b2*beta^2*ds/2)*dq;
    
    % Non-circulatory - AOA contribution
    Kprime_a(i) = Kprime_a(i-1)*exp(-dt/T_alpha)+(dat(i)-dat(i-1))*exp(-dt/(2*T_alpha));
    
    % Non-circulatory - Pitch rate contribution
    Kprime_q(i) = Kprime_q(i-1)*exp(-dt/T_q)+(dqt(i)-dqt(i-1))*exp(-dt/(2*T_q));
    
end

% Circulatory normal force
alpha_E= alpha-alpha0-X1-X2;
Cnc = C_Nalphac*alpha_E;

% Non-circulatory normal force
Cni_a = 4*T_alpha/M*(dat-Kprime_a);
Cni_q = -T_q/M*(dqt-Kprime_q);
Cni = Cni_a+Cni_q;

% Total normal force under attached conditions
Cnp = Cnc+Cni;


% AERODYNAMIC MOMENT COEFFICIENT
for i = 2:N
    
    ds = s_span(i)-s_span(i-1);
    dq = q(i)-q(i-1);
    dt = t_span(i)-t_span(i-1);
    
    % Circulatory - Pitch rate contribution
    Kprimeprime_qM(i) = Kprimeprime_qM(i-1)*exp(-b5*beta^2*ds)+A5*dq*exp(-b5*beta^2*ds/2);
    
    % Non-circulatory - Pitch rate contribution
    Kprime_qM(i) = Kprime_qM(i-1)*exp(-dt/(K_qM^2*T_I))+(dqt(i)-dqt(i-1))*exp(-dt/(2*K_qM^2*T_I));
    
end

% Circulatory moment
Cmc_a = -C_Nalphac*(1-A1*exp(-b1*beta^2*s_span)-A2*exp(-b2*beta^2*s_span))*(x_ac-0.25);
Cmc_q = -C_Nalpha/(16*beta)*(q-Kprimeprime_qM)*c/V;
Cmc = Cmc_a+Cmc_q;

% Non-circulatory moment
Cmi_a = -Cni_a/4;
Cmi_q = -7*K_qM^2*T_I/(12*M)*(dqt-Kprime_qM);
Cmi = Cmi_a+Cmi_q;

% Total aerodynamic moment coefficient
Cmp = Cmc+Cmi;


% CHORD FORCE
Ccp = Cnc.*tan(alpha_E+alpha0);

end