%% Simulation of the Beddoes-Leishman State Space Model (NACA0012)
% University of California, Irvine - Fall 2021
% Laura Pla Olea - lplaolea@uci.edu

clear; clc; close all;

%% Input data

% Geometry
airfoil = ('NACA0012');
c = 1; % Airfoil chord [m]
e = -1/2*c/2; % Pitching axis (measured from the midchord) [m]
M = 0.3; % Mach number
k = 0.1; % Reduced frequency
n_cycles = 10; % Number of cycles that are computed
n_t = 1e2; % Number of time steps per cycle

% Motion
alphabase = deg2rad(5); % Initial angle [rad]
A_alpha = deg2rad(10); % Pitching amplitude [rad]
H = 0; % Plunging amplitude (h/c)
phi = deg2rad(0); % Phase between the pitching and plunging motions [rad]

% Environment
gamma = 1.4; % Heat capacity ratio / Adiabatic index
R = 287.058; % Specific gas constant for dry air [J/kg·K]
Ta = 298; % Ambient temperature [K]

% Numerical values
[eta, C_Nalpha, alpha0, Cd0, Cm0, alpha1, dalpha1, S1, S2, K0, K1, K2, T_P, T_f, C_N1, T_v, T_vl, D_f] = input_NACA0012(M);
x_ac = 0.25-K0; % Aerodynamic center normalized by the chord

%% Preliminary calculations

a = sqrt(gamma*R*Ta); % Speed of sound [m/s]
V = M*a; % Free-stream velocity [m/s]

b = c/2; % Half-chord [m]
w = V*k/b; % Angular velocity [rad/s]
T = 2*pi/w; % Period [s]


%% Kinematics

N = n_cycles*n_t;
t = linspace(0,n_cycles*T,N); % Time vector [s]
s = V*t/b; % Non-dimensional time vector

% Pitching motion [rad]
alpha = alphabase-A_alpha*cos(w*t);
dalpha = w*A_alpha*sin(w*t); % First derivative [rad/s]

% Plunging motion [m]
dh = w*H*b*sin(w*t+phi);

% Effective angle of attack [rad]
alpha_eff = alpha+atan(dh./V);

%% Attached flow

q = dalpha*c/V; % Non-dimensional pitch rate

[Cnp, Cmp, Ccp, alpha_E, Cni] = BL_attached(t,alpha_eff,q,V,M,c,C_Nalpha,alpha0,Cm0,x_ac);

%% Stall onset

Cnprime = BL_stallonset(s,T_P,C_Nalpha*alpha_E,Cnp);

%% Trailing edge separation

[Cnf, Cmf, Ccf, fprimeprime, fprime] = BL_TEseparation(s,Cnprime,C_Nalpha,C_N1,alpha0,alpha,alpha_E,alpha1,dalpha1,S1,S2,T_f,T_vl,K0,K1,K2,eta,D_f);

%% Modeling of dynamic stall

[Cnv, Cmv] = BL_dynamicstall(s,fprimeprime,alpha,alpha_E,Cnprime,C_Nalpha,T_v,T_vl,C_N1);

%% Final results

% Total normal force
CN = Cnf+Cnv+Cni;
CM = Cmf+Cmv+Cmp-K0*C_Nalpha*alpha_E;

figure(1);
plot(rad2deg(alpha_eff(N-n_t:N)),Cnp(N-n_t:N),rad2deg(alpha_eff(N-n_t:N)),Cnf(N-n_t:N),rad2deg(alpha_eff(N-n_t:N)),Cnv(N-n_t:N));
xlabel('\alpha [º]'); legend('C_{N}^{p}','C_{N}^{f}','C_{N}^{v}','Location','Northwest');
% axis ([-10 20 -1 2]);
grid on
title([airfoil,', k=',num2str(k),', $\alpha=',num2str(rad2deg(alphabase)),'^{o}+',num2str(rad2deg(A_alpha)),'^{o}sin(\omega t)$'],'interpreter','latex');

% Lift & Drag coefficients
CL = CN.*cos(alpha_eff)+Ccf.*sin(alpha_eff);
CD = CN.*sin(alpha_eff)-Ccf.*cos(alpha_eff)+Cd0;

% Comparison with experimental results
[alpha_el,CN_e,alpha_ed,CD_e,alpha_em,CM_e] = experimental_results(airfoil,k,alphabase,A_alpha);

figure(2);
plot(rad2deg(alpha_eff(N-n_t:N)),CN(N-n_t:N),'r')
xlabel('\alpha [º]'); ylabel('C_{N}');
grid on
hold on; plot (alpha_el,CN_e,'--')
hold on; plot (alpha_el,CN_e,'o')
legend('Model','Experimental','Location','best');
title([airfoil,', k=',num2str(k),', $\alpha=',num2str(rad2deg(alphabase)),'^{o}+',num2str(rad2deg(A_alpha)),'^{o}sin(\omega t)$'],'interpreter','latex');

figure(3);
plot(rad2deg(alpha_eff(N-n_t:N)),CD(N-n_t:N),'r')
grid on
hold on; plot(alpha_ed,CD_e,'--')
hold on; plot (alpha_ed,CD_e,'o')
legend('Model','Experimental','Location','best')
xlabel('\alpha [º]');
ylabel('C_{D}');
title([airfoil,', k=',num2str(k),', $\alpha=',num2str(rad2deg(alphabase)),'^{o}+',num2str(rad2deg(A_alpha)),'^{o}sin(\omega t)$'],'interpreter','latex');

figure(4);
plot(rad2deg(alpha_eff(N-n_t:N)),CM(N-n_t:N),'r')
grid on
hold on; plot(alpha_em,CM_e,'--')
hold on; plot (alpha_em,CM_e,'o')
legend('Model','Experimental','Location','best')
xlabel('\alpha [º]');
ylabel('C_{M}'); xlim([-5 15]); ylim([-.2 .1]);
title([airfoil,', k=',num2str(k),', $\alpha=',num2str(rad2deg(alphabase)),'^{o}+',num2str(rad2deg(A_alpha)),'^{o}sin(\omega t)$'],'interpreter','latex');