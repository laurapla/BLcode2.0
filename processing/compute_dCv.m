% Simulation of the Beddoes-Leishman State Space Model (NACA0012) to find
% the value of dCv as a function of time and alphabase
% University of California, Irvine - Fall 2021
% Laura Pla Olea - lplaolea@uci.edu

clear; clc; close all
addpath(genpath('../'))

%% Input data

% Geometry
airfoil = ('NACA0012');
c = 1; % Airfoil chord [m]
e = -1/2*1/2; % Pitching axis (measured from the midchord) normalized by the chord
M = 0.3; % Mach number
k = 0.3; % Reduced frequency
n_cycles = 10; % Number of cycles that are computed
n_t = 1e2; % Number of time steps per cycle

% Motion
A_alpha = deg2rad(0); % Pitching amplitude [rad]
H = 0.55; % Plunging amplitude (h/c)
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
dalpha = w*A_alpha*sin(w*t); % First derivative [rad/s]

% Plunging motion [m]
dh = w*H*b*sin(w*t+phi);

%% Beddoes-Leishman Model

alpha_array = linspace(0,pi/4,200);
dCv_array = zeros(n_t+1,200);

for j = 1:200
    
    % Pitching motion [rad]
    alphabase = alpha_array(j);
    alpha = alphabase-A_alpha*cos(w*t);
    
    % Effective angle of attack [rad]
    alpha_eff = alpha+atan(dh/V);
    
    % Attached flow
    
    q = dalpha*c/V; % Non-dimensional pitch rate
    
    [Cnp, Cmp, Ccp, alpha_E, Cni] = BL_attached(t,alpha_eff,q,V,M,c,e,C_Nalpha,alpha0,Cm0,x_ac);
    
    % Stall onset
    
    Cnprime = BL_stallonset(s,T_P,C_Nalpha*alpha_E,Cnp);
    
    % Trailing edge separation
    
    [Cnf, Cmf, Ccf, fprimeprime, fprime] = BL_TEseparation(s,Cnprime,C_Nalpha,C_N1,alpha0,alpha_eff,alpha_E,alpha1,dalpha1,S1,S2,T_f,T_vl,K0,K1,K2,eta,D_f);
    
    % Modeling of dynamic stall
    
    % Time rate of change of circulation
    tauv = vortex_time(Cnprime,s,C_N1);
    
    Cnc = C_Nalpha*alpha_E; % Circulatory normal force coefficient
    Kn = (1+sqrt(fprimeprime)).^2/4;
    
    N = size(s,2);
    
    Cv = zeros(1,N);
    Cnv = zeros(1,N);
    
    for i = 2:N
        
        ds_span = s(i)-s(i-1);
        Sa = alpha_eff(i)-alpha_eff(i-1);
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
        
    end
    
    dCv_array(:,j) = gradient(Cv(N-n_t:N))./gradient(t(N-n_t:N));
    
end

filename = ['dCv_H' num2str(100*H)];
save(filename,'alpha_array','dCv_array','t');