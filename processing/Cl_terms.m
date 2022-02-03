%% Separately computing the terms of the lift coefficient averaged over a cycle
% University of California, Irvine - Fall 2021
% Laura Pla Olea - lplaolea@uci.edu

% clear; clc; close all;

%% Loading files

addpath(genpath('../'),genpath('/data'),genpath('../experimental_data'))
load('dCv_aA.mat');
load('CN_NACA0012_static.mat');

%% Input data

% Geometry
airfoil = ('NACA0012');
M = 0.3; % Mach number
k = 0.1; % Reduced frequency

% Motion
n_A = 16;
A_alpha = linspace(0,15,n_A); % Pitching amplitude [rad]
H = 0; % Plunging amplitude (h/c)
phi = deg2rad(0); % Phase between the pitching and plunging motions [rad]

% Numerical values
[eta, C_Nalpha, ~, ~, ~, alpha1, ~, S1, S2, ~, ~, ~, ~, ~, ~, ~, ~, ~] = input_NACA0012(M);

% Angle parameters
N = 200;
alpha = deg2rad(linspace(0,45,N)).';

%% Averaged Cl

% State x10 (separation point)
x10 = separation_point(alpha-H*k*sin(phi),alpha1,S1,S2).';

% State x11 (dCv)
[oldcols,oldrows] = meshgrid(A_array,alpha_array);
[newcols, newrows] = meshgrid(A_alpha,alpha);
x11 = interp2(oldcols,oldrows,dCv_alphaA,newcols,newrows);

%% Averaged Cl

t1 = zeros(N,n_A);
t2 = zeros(N,n_A);
t3 = zeros(N,n_A);
t4 = zeros(N,n_A);
t5 = zeros(N,n_A);
t6 = zeros(N,n_A);
t7 = zeros(N,n_A);
t8 = zeros(N,n_A);

for j = 1:n_A
    
    % Averaged lift coefficient
    t1(:,j) = x11(:,j).*cos(alpha(i));
    t2(:,j) = C_Nalpha/2*(alpha-H*k*sin(phi)).*(1+sqrt(x10)).^2/2.*cos(alpha);
    t3(:,j) = C_Nalpha/2*(alpha-H*k*sin(phi))*2*eta.*(alpha-H*k*sin(phi)).*sqrt(x10).*sin(alpha);
    t4(:,j) = -deg2rad(A_alpha(j))^2/4*x11(:,j).*cos(alpha);
    t5(:,j) = -deg2rad(A_alpha(j))^2/4*sin(alpha)/(8*M);
    t6(:,j) = -deg2rad(A_alpha(j))^2/4*C_Nalpha*alpha/2.*(1+sqrt(x10)).^2/2.*cos(alpha);
    t7(:,j) = -deg2rad(A_alpha(j))^2/4*eta*C_Nalpha*alpha.^2.*sqrt(x10).*sin(alpha);
    t8(:,j) = deg2rad(A_alpha(j))*H*k/M*sin(alpha)*sin(phi);
    
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
        Legend{iter} = strcat('A_{\alpha}=', num2str((iter-2)), 'º');
    end
end
legend(Legend,'Location','bestoutside')
grid on; grid minor

%% Terms

index = 6;
figure;
plot(rad2deg(alpha),t1(:,index),rad2deg(alpha),t2(:,index),rad2deg(alpha),t3(:,index),rad2deg(alpha),t4(:,index),rad2deg(alpha),t5(:,index),rad2deg(alpha),t6(:,index),rad2deg(alpha),t7(:,index),rad2deg(alpha),t8(:,index),rad2deg(alpha),Cl(:,index));
legend('$x_{11}^{*}\cos\alpha^{*}$','$C_{N_{\alpha}}(\alpha^{*}-Hbk\sin\phi)(1+\sqrt{x_{10}^{*}})^2\cos\alpha^{*}/4$'...
    ,'$\eta C_{N_{\alpha}}(\alpha^{*}-Hbk\sin\phi)^2\sqrt{x_{10}^{*}}\sin\alpha^{*}$'...
    ,'$-A_{\alpha}^{2}x_{11}^{*}\cos\alpha^{*}/4$','$-A_{\alpha}^{2}\sin\alpha^{*}/32M$'...
    ,'$-A_{\alpha}^{2}C_{N_{\alpha}}\alpha^{*}(1+\sqrt{x_{10}^{*}})^2\cos\alpha^{*}/16$'...
    ,'$-A_{\alpha}^{2}\eta C_{N_{\alpha}}\alpha^{*2}\sqrt{x_{10}^{*}}\sin\alpha^{*}/4$'...
    ,'$2A_{\alpha}Hk\sin\alpha^{*}\sin\phi/M$','Cl','interpreter','latex')