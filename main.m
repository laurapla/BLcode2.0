%% Simulation of the Beddoes-Leishman State Space Model (NACA0012)
% University of California, Irvine - Fall 2021
% Laura Pla Olea - lplaolea@uci.edu

clear; clc; close all;

%% Input data

% Geometry
airfoil = ('NACA0012');
c = 1; % Airfoil chord [m]
e = -1/2*1/2; % Pitching axis (measured from the midchord) normalized by the chord
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
R = 287.058; % Specific gas constant for dry air [J/kg�K]
Ta = 298; % Ambient temperature [K]

% Numerical values
[eta, C_Nalpha, alpha0, Cd0, Cm0, alpha1, dalpha1, S1, S2, K0, K1, K2, T_P, T_f, C_N1, T_v, T_vl, D_f] = input_NACA0012(M);
% C_Nalpha = 5.9303;
% alpha1 = deg2rad(15.0000);
% S1 = deg2rad(5.2849);
% S2 = deg2rad(1.0337);
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

[Cnp, Cmp, Ccp, alpha_E, Cni] = BL_attached(t,alpha_eff,q,V,M,c,e,C_Nalpha,alpha0,Cm0,x_ac);

%% Stall onset

Cnprime = BL_stallonset(s,T_P,Cnp);

%% Trailing edge separation

[Cnf, Cmf, Ccf, fprimeprime] = BL_TEseparation(s,Cnprime,Cni,C_Nalpha,C_N1,alpha,alpha_E,alpha0,alpha1,dalpha1,S1,S2,T_f,T_vl,K0,K1,K2,eta,D_f);

%% Modeling of dynamic stall

[Cnv, Cmv] = BL_dynamicstall(s,fprimeprime,alpha,alpha_E,Cnprime,C_Nalpha,T_v,T_vl,C_N1);

%% Final results

% Total normal force
CN = Cnf+Cnv;
CM = Cmf+Cmv+Cmp+C_Nalpha*alpha_E*(x_ac-0.25);

lwidth = 1.7;
lbl_font = 14;
lgd_font = 14;
ax_font = 12;

figure(1);
plot(rad2deg(alpha_eff(N-n_t:N)),Cnp(N-n_t:N),rad2deg(alpha_eff(N-n_t:N)),Cnf(N-n_t:N),rad2deg(alpha_eff(N-n_t:N)),Cnv(N-n_t:N),'LineWidth',lwidth);
xlabel('$\alpha, ^{\circ}$','Interpreter','latex','FontSize',lbl_font);
legend('$C_{N}^{p}$','$C_{N}^{f}$','$C_{N}^{v}$','Location','Northwest','Interpreter','latex','FontSize',lgd_font);
grid on
title([airfoil,', $k=',num2str(k),'$, $\alpha=',num2str(rad2deg(alphabase)),'^{o}+',num2str(rad2deg(A_alpha)),'^{o}sin(\omega t)$'],'interpreter','latex');

% Lift & Drag coefficients
CL = CN.*cos(alpha_eff)+Ccf.*sin(alpha_eff);
CD = CN.*sin(alpha_eff)-Ccf.*cos(alpha_eff)+Cd0;

% Comparison with experimental results
[alpha_el,CL_e,alpha_ed,CD_e,alpha_em,CM_e] = experimental_results(airfoil,k,M,alphabase,A_alpha,H,phi);

figure(2);
plot(rad2deg(alpha_eff(N-n_t:N)),CN(N-n_t:N),'r','LineWidth',lwidth)
xlabel('$\alpha, ^{\circ}$','Interpreter','latex','FontSize',lbl_font);
ylabel('$C_{n}$','Interpreter','latex','FontSize',lbl_font);
grid on
hold on; plot(alpha_el,CL_e,'--o','LineWidth',lwidth,'Color',[0.9290 0.6940 0.1250])
legend('Model','Experimental','Location','best','FontSize',lgd_font);
title([airfoil,', k=',num2str(k),', $\alpha=',num2str(rad2deg(alphabase)),'^{o}+',num2str(rad2deg(A_alpha)),'^{o}sin(\omega t)$'],'interpreter','latex');
ax = gca;
ax.FontSize = ax_font;

figure(3);
plot(rad2deg(alpha_eff(N-n_t:N)),CD(N-n_t:N),'r','LineWidth',lwidth)
grid on
hold on; plot(alpha_ed,CD_e,'--o','LineWidth',lwidth,'Color',[0.9290 0.6940 0.1250])
legend('Model','Experimental','Location','best','FontSize',lgd_font)
xlabel('$\alpha, ^{\circ}$','Interpreter','latex','FontSize',lbl_font);
ylabel('$C_{D}$','Interpreter','latex','FontSize',lbl_font);
title([airfoil,', k=',num2str(k),', $\alpha=',num2str(rad2deg(alphabase)),'^{o}+',num2str(rad2deg(A_alpha)),'^{o}sin(\omega t)$'],'interpreter','latex');
ax = gca; % current axes
ax.FontSize = ax_font;

figure(4);
plot(rad2deg(alpha_eff(N-n_t:N)),CM(N-n_t:N),'r','LineWidth',lwidth)
grid on
hold on; plot(alpha_em,CM_e,'--o','LineWidth',lwidth,'Color',[0.9290 0.6940 0.1250])
legend('Model','Experimental','Location','best','FontSize',lgd_font)
xlabel('$\alpha, ^{\circ}$','Interpreter','latex','FontSize',lbl_font);
ylabel('$C_{M}$','Interpreter','latex','FontSize',lbl_font);
title([airfoil,', k=',num2str(k),', $\alpha=',num2str(rad2deg(alphabase)),'^{o}+',num2str(rad2deg(A_alpha)),'^{o}sin(\omega t)$'],'interpreter','latex');
ax = gca; % current axes
ax.FontSize = ax_font;