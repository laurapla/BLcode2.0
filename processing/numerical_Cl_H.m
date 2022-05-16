% Simulation of the Beddoes-Leishman State Space Model (NACA0012) to find
% the value of dCv as a function of time and alphabase
% University of California, Irvine - Winter 2022
% Laura Pla Olea - lplaolea@uci.edu

clear; clc; close all
addpath(genpath('../'))

%% Input data

% Geometry
airfoil = ('NACA0012');
c = 1; % Airfoil chord [m]
e = -1/2*1/2; % Pitching axis (measured from the midchord) normalized by the chord
M = 0.3; % Mach number
k = 0.5; % Reduced frequency
n_cycles = 10; % Number of cycles that are computed
n_t = 1e2; % Number of time steps per cycle

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


%% Kinematics + Beddoes-Leishman model

N = n_cycles*n_t;
t = linspace(0,n_cycles*T,N); % Time vector [s]
s = V*t/b; % Non-dimensional time vector

n_H = 12;
n_alpha = 20;

H_array = linspace(0,0.55,n_H);
alpha_array = linspace(0,pi/4,n_alpha); % Mean angle of attack

avg_Cl = zeros(n_alpha,n_H);
avg_Cd = zeros(n_alpha,n_H);

for index = 1:n_H
    
    % Motion
    A_alpha = deg2rad(0); % Pitching amplitude [rad]
    H = H_array(index); % Plunging amplitude (h/c)
    phi = deg2rad(0); % Phase between the pitching and plunging motions [rad]
    
    % Pitching motion [rad]
    dalpha = w*A_alpha*sin(w*t); % First derivative [rad/s]
    
    % Plunging motion [m]
    dh = w*H*b*sin(w*t+phi);
    
    % Beddoes-Leishman Model
        
    for j = 1:n_alpha
        
        % Pitching motion [rad]
        alphabase = alpha_array(j);
        alpha = alphabase-A_alpha*cos(w*t);
        
        % Effective angle of attack [rad]
        alpha_eff = alpha+atan(dh/V);
        
        
        % Attached flow
        
        q = dalpha*c/V; % Non-dimensional pitch rate
        
        [Cnp, ~, ~, alpha_E, Cni] = BL_attached(t,alpha_eff,q,V,M,c,e,C_Nalpha,alpha0,Cm0,x_ac);
        
        
        % Stall onset
        
        Cnprime = BL_stallonset(s,T_P,C_Nalpha*alpha_E,Cnp);
        
        
        % Trailing edge separation
        
        [Cnf, ~, Ccf, fprimeprime, ~] = BL_TEseparation(s,Cnprime,C_Nalpha,C_N1,alpha0,alpha_eff,alpha_E,alpha1,dalpha1,S1,S2,T_f,T_vl,K0,K1,K2,eta,D_f);
        
        
        % Modeling of dynamic stall
        
        [Cnv, Cmv] = BL_dynamicstall(s,fprimeprime,alpha,alpha_E,Cnprime,C_Nalpha,T_v,T_vl,C_N1);


        % Total force
        
        CN = Cnf+Cnv+Cni;
        CL = CN.*cos(alpha_eff)+Ccf.*sin(alpha_eff);
        CD = CN.*sin(alpha_eff)-Ccf.*cos(alpha_eff)+Cd0;


        % Average values

        avg_Cl(j,index) = mean(CL(N-n_t:N));
        avg_Cd(j,index) = mean(CD(N-n_t:N));
        
    end
    
end

filename = ['Cl_Cd_H_k',num2str(100*k)];
save(filename,'H_array','alpha_array','avg_Cl','avg_Cd');

%% Plot average CL vs. alpha (for different H)

width = 1.7;
font_lgd = 10;
font_labels = 14;

symbols = {'-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.'};
Okabe_Ito = [0.902 0.624 0; 0.337 0.737 0.914; 0 0.62 0.451;
0.941 0.894 0.259; 0 0.447 0.698; 0.835 0.369 0; 0.8 0.475 0.655];

% Plotting the results
interval = 1;
Legend_H = cell(n_H,1);

figure(1);
colororder(Okabe_Ito)
for j = 1:n_H
    if rem(100*H_array(j),interval)<1e-10
        hold on;
        if j==1
            plot(rad2deg(alpha_array),avg_Cl(:,j),'k','LineWidth',width)
        else
            plot(rad2deg(alpha_array),avg_Cl(:,j),symbols{j},'LineWidth',width)
        end
        Legend_H{j} = strcat('$H=$', num2str(H_array(j)));
    end
end

legend(Legend_H,'Location','best','interpreter','latex','FontSize',font_lgd)
xlabel('$\alpha^{*}, ^{\circ}$','interpreter','latex','FontSize',font_labels);
ylabel('$\overline{C}_{L}$','interpreter','latex','FontSize',font_labels);
grid on;

%% Plot average CD vs. alpha (for different H)

figure(2);
colororder(Okabe_Ito)
for j = 1:n_H
    if rem(100*H_array(j),interval)<1e-10
        hold on;
        if j==1
            plot(rad2deg(alpha_array),avg_Cd(:,j),'k','LineWidth',width)
        else
            plot(rad2deg(alpha_array),avg_Cd(:,j),symbols{j},'LineWidth',width)
        end
    end
end

legend(Legend_H,'Location','best','interpreter','latex','FontSize',font_lgd)
xlabel('$\alpha^{*}, ^{\circ}$','interpreter','latex','FontSize',font_labels);
ylabel('$\overline{C}_{D}$','interpreter','latex','FontSize',font_labels);
grid on;