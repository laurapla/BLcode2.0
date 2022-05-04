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
k = 0.5; % Reduced frequency
n_cycles = 10; % Number of cycles that are computed
n_t = 1e2; % Number of time steps per cycle

% Environment
gamma = 1.4; % Heat capacity ratio / Adiabatic index
R = 287.058; % Specific gas constant for dry air [J/kg�K]
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

n_A = 11;
n_alpha = 20;
A_array = linspace(0,10,n_A);
dCv_matrix = zeros(n_t+1,n_alpha,n_A);

for index = 1:n_A
    
    % Motion
    A_alpha = deg2rad(A_array(index)); % Pitching amplitude [rad]
    H = 0; % Plunging amplitude (h/c)
    phi = deg2rad(0); % Phase between the pitching and plunging motions [rad]
    
    % Pitching motion [rad]
    dalpha = w*A_alpha*sin(w*t); % First derivative [rad/s]
    
    % Plunging motion [m]
    dh = w*H*b*sin(w*t+phi);
    
    % Beddoes-Leishman Model
    
    % Mean angle of attack
    alpha_array = linspace(0,pi/4,n_alpha);
    
    for j = 1:n_alpha
        
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
        Ds = zeros(1,N);
        
        for j_index = 2:N
            
            ds_span = s(j_index)-s(j_index-1);
            Sa = alpha_eff(j_index)-alpha_eff(j_index-1);
            df = fprimeprime(j_index)-fprimeprime(j_index-1);
            sigma = sigma2(tauv(j_index), T_vl, Sa, df);
            
            if tauv(j_index)>0 && tauv(j_index)<=T_vl
                Ds(j_index) = 1;
            elseif abs(Cnprime(j_index))<C_N1 && df<0
                Ds(j_index) = 1;
            elseif Sa>0 && df>0
                Ds(j_index) = 1;
            else
                Ds(j_index) = 0;
            end
            
            Cv(j_index) = Ds(j_index)*Cnc(j_index)*(1-Kn(j_index));
            
            Ds2 = Ds;
            if Ds(j_index)==1 && Ds(j_index-1)==0
                Ds2(j_index)=0;
            elseif j_index<N
                if Ds(j_index)==1 && Ds(j_index+1)==0
                    Ds2(j_index)=0;
                end
            end
            
        end
        
        dCv_matrix(:,j,index) = Ds2(N-n_t:N).*gradient(Cv(N-n_t:N))./gradient(t(N-n_t:N));
        
    end
    
end

filename = ['dCv_A_k',num2str(100*k)];
save(filename,'A_array','dCv_matrix','t');

%% Processing

% Storing the results
dCv = zeros(n_alpha,n_A);
interval = 2;

for j = 1:n_A
    
    dCv_array = dCv_matrix(:,:,j);
    
    % Storing the results
    dCv(:,j) = mean(dCv_array);
    dCv(isinf(dCv)|isnan(dCv)) = 0;
    
    if rem(round(A_array(j)),interval)==0
        hold on;
        plot(rad2deg(alpha_array),dCv(:,j))
    end
    
end

%% Plot parameters

Legend = cell(floor(n_A/interval),1);
for iter=1:n_A
    if rem(A_array(iter),interval)==0
        Legend{(ceil(iter/interval))} = strcat('A_{\alpha}=', num2str((A_array(iter))), '�');
    end
end
legend(Legend,'Location','northwest')

xlabel('\alpha [�]'); ylabel('$\dot{C}_{v}^{*}$','interpreter','latex');

% %% Plot as a function of the parameter
% 
% interval_angle = 1;
% fffit = zeros(3,floor(n_alpha/interval_angle));
% 
% figure(2);
% for j = 1:n_alpha
%     f = fit(A_array.',dCv(j,:).','poly2');
%     fffit(:,j) = coeffvalues(f);
%     if rem(j-1,interval_angle)==0
%         hold on;
%         plot(A_array,dCv(j,:));
%     end
% end
% 
% Legend2 = cell(floor(n_alpha/interval_angle),1);
% for iter = 1:n_alpha
%     if rem(iter-1,interval_angle)==0
%         Legend2{(iter-1)/interval_angle+1} = strcat('\alpha^{*}=', num2str(round(rad2deg(alpha_array(iter)))), '�');
%     end
% end
% legend(Legend2,'Location','bestoutside','NumColumns',2)
% 
% xlabel('A_{\alpha} [�]');
% 
% ylabel('$\dot{C}_{v}^{*}$','interpreter','latex');