%% Separately computing the terms of the lift coefficient averaged over a cycle
% University of California, Irvine - Fall 2021
% Laura Pla Olea - lplaolea@uci.edu

clear; clc; close all;

%% Input data

% Geometry
airfoil = ('NACA0012');
M = 0.3; % Mach number
k = 0.1; % Reduced frequency

%% Loading files

addpath(genpath('../'),genpath('/data'),genpath('../experimental_data'))
load(strcat('dCvM', num2str(M), '.mat'));
load('CN_NACA0012_static.mat');

%% Pre-calculations

n_A_array = length(A_array);
n_k_array = length(k_array);
n_alpha_array = length(alpha_array);

% Motion
A_alpha = linspace(0,10,n_A_array); % Pitching amplitude [rad]
H = 0; % Plunging amplitude (h/c)
phi = deg2rad(0); % Phase between the pitching and plunging motions [rad]

% Numerical values
[eta, C_Nalpha, ~, ~, ~, alpha1, ~, S1, S2, ~, ~, ~, ~, ~, ~, ~, ~, ~] = input_NACA0012(M);

% Angle parameters
N = 200;
alpha = deg2rad(linspace(0,45,N)).';

%% dCv_matrix processing

k_index = 1;
for i = 2:n_k_array
    if abs(k-k_array(i))<abs(k-k_array(i-1))
        k_index = i;
    end
end

%% States x10 and x11

% State x10 (separation point)
x10 = separation_point(alpha-H*k*sin(phi),alpha1,S1,S2).';

% State x11 (dCv)
[oldcols,oldrows] = meshgrid(A_array,alpha_array);
[newcols, newrows] = meshgrid(A_alpha,alpha);
x11 = interp2(oldcols,oldrows,dCv_A(:,:,k_index),newcols,newrows);

%% Average Cl

t1 = zeros(N,n_A_array);
t2 = zeros(N,n_A_array);
t3 = zeros(N,n_A_array);
t4 = zeros(N,n_A_array);
t5 = zeros(N,n_A_array);
t6 = zeros(N,n_A_array);
t7 = zeros(N,n_A_array);
t8 = zeros(N,n_A_array);

for j = 1:n_A_array
    
    % Averaged lift coefficient
    t1(:,j) = x11(:,j).*cos(alpha);
    t2(:,j) = C_Nalpha/2*(alpha-H*k*sin(phi)).*(1+sqrt(x10)).^2/2.*cos(alpha);
    t3(:,j) = C_Nalpha/2*(alpha-H*k*sin(phi))*2*eta.*(alpha-H*k*sin(phi)).*sqrt(x10).*sin(alpha);
    t4(:,j) = -deg2rad(A_alpha(j))^2/4*x11(:,j).*cos(alpha);
    t5(:,j) = -deg2rad(A_alpha(j))^2/4*8*sin(alpha)/M;
    t6(:,j) = -deg2rad(A_alpha(j))^2/4*C_Nalpha*alpha/2.*(1+sqrt(x10)).^2/2.*cos(alpha);
    t7(:,j) = -deg2rad(A_alpha(j))^2/4*eta*C_Nalpha*alpha.^2.*sqrt(x10).*sin(alpha);
    t8(:,j) = 2*deg2rad(A_alpha(j))*H*k/M*sin(alpha)*sin(phi);
    
end

Cl = t1+t2+t3+t4+t5+t6+t7+t8;

%% Static value

alpha_s = CN_NACA0012_static(:,1);
Cl_s = CN_NACA0012_static(:,2);

figure;
plot(alpha_s,Cl_s,'--o');
for i = 1:8
    hold on;
    plot(rad2deg(alpha),Cl(:,i));
end

xlabel('\alpha^{*}'); ylabel('$\overline{C}_{L}$','interpreter','latex');

Legend = cell(8,1);
for iter = 1:8+1
    if iter==1
        Legend{iter} = strcat('Steady Cl');
    else
        Legend{iter} = strcat('A_{\alpha}=', num2str((iter-2)), '�');
    end
end
legend(Legend,'Location','bestoutside')
grid on; grid minor

%% Terms

index = 6;
figure;
plot(rad2deg(alpha),t1(:,index),rad2deg(alpha),t2(:,index),rad2deg(alpha),t3(:,index),rad2deg(alpha),t4(:,index),rad2deg(alpha),t5(:,index),rad2deg(alpha),t6(:,index),rad2deg(alpha),t7(:,index),rad2deg(alpha),t8(:,index));
legend('$x_{11}^{*}\cos\alpha^{*}$','$C_{N_{\alpha}}(\alpha^{*}-Hbk\sin\phi)(1+\sqrt{x_{10}^{*}})^2\cos\alpha^{*}/4$'...
    ,'$\eta C_{N_{\alpha}}(\alpha^{*}-Hbk\sin\phi)^2\sqrt{x_{10}^{*}}\sin\alpha^{*}$'...
    ,'$-A_{\alpha}^{2}x_{11}^{*}\cos\alpha^{*}/4$','$-2A_{\alpha}^{2}\sin\alpha^{*}/M$'...
    ,'$-A_{\alpha}^{2}C_{N_{\alpha}}\alpha^{*}(1+\sqrt{x_{10}^{*}})^2\cos\alpha^{*}/16$'...
    ,'$-A_{\alpha}^{2}\eta C_{N_{\alpha}}\alpha^{*2}\sqrt{x_{10}^{*}}\sin\alpha^{*}/4$'...
    ,'$2A_{\alpha}Hk\sin\alpha^{*}\sin\phi/M$','interpreter','latex')