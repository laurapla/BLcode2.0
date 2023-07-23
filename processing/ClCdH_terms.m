%% Separately computing the terms of the lift coefficient averaged over a cycle
% University of California, Irvine - Fall 2022
% Laura Pla Olea - lplaolea@uci.edu

clear; clc; close all;

%% Input data

% Geometry
airfoil = ('NACA0012');
M = 0.7; % Mach number
k = 0.5; % Reduced frequency
H_display = 0.2; % Amplitude to look at

%% Loading files

addpath(genpath('/dCv_data'))
load(strcat('dCvM', num2str(M), '.mat'));

%% Pre-calculations

n_H_array = length(H_array);
n_k_array = length(k_array);
n_alpha_array = length(alpha_array);

% Motion
A_alpha = 0; % Pitching amplitude [rad]
H = linspace(0,0.55,n_H_array); % Plunging amplitude (normalized by the half-chord)
phi = deg2rad(0); % Phase between the pitching and plunging motions [rad]

% Numerical values
[eta, C_Nalpha, ~, ~, ~, alpha1, ~, S1, S2, ~, ~, ~, ~, ~, ~, ~, ~, ~] = input_NACA0012(M);
C_Nalpha = 5.9303;
alpha1 = deg2rad(15);
S1 = deg2rad(5.284944901391905);
S2 = deg2rad(1.033663269609409);

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
x10 = separation_point(alpha,alpha1,S1,S2).';

% State x11 (dCv)
[oldcols,oldrows] = meshgrid(H_array,alpha_array);
[newcols, newrows] = meshgrid(H,alpha);
x11 = interp2(oldcols,oldrows,dCv_H(:,:,k_index),newcols,newrows);

%% Average Cl

n_terms = 3;
term_Cl = zeros(N,n_H_array,n_terms);
term_Cd = zeros(N,n_H_array,n_terms);

for j = 1:n_H_array
    
    % Average lift coefficient
    term_Cl(:,j,1) = x11(:,j).*cos(alpha);
    term_Cl(:,j,2) = C_Nalpha/2*(alpha-H(j)*k*sin(phi)).*(1+sqrt(x10)).^2/2.*cos(alpha);
    term_Cl(:,j,3) = C_Nalpha/2*(alpha-H(j)*k*sin(phi))*2*eta.*(alpha-H(j)*k*sin(phi)).*sqrt(x10).*sin(alpha);

    % Average drag coefficient
    term_Cd(:,j,1) = x11(:,j).*sin(alpha);
    term_Cd(:,j,2) = C_Nalpha/2*(alpha-H(j)*k*sin(phi)).*(1+sqrt(x10)).^2/2.*sin(alpha);
    term_Cd(:,j,3) = -C_Nalpha/2*(alpha-H(j)*k*sin(phi))*2*eta.*(alpha-H(j)*k*sin(phi)).*sqrt(x10).*cos(alpha);
    
end

Cl = zeros(N,n_H_array);
Cd = zeros(N,n_H_array);
for i = 1:n_terms
    Cl = Cl+term_Cl(:,:,i);
    Cd = Cd+term_Cd(:,:,i);
end

%% Total value

line_width = 1.7;
font_lgd = 10;
font_labels = 14;
inter = 2;

symbols = {'-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.'};
Okabe_Ito = [0.902 0.624 0; 0.337 0.737 0.914; 0 0.62 0.451;
    0.941 0.894 0.259; 0 0.447 0.698; 0.835 0.369 0; 0.8 0.475 0.655];

% Cl

figure;
colororder(Okabe_Ito)
for i = 1:n_H_array
    if i==1
        plot(rad2deg(alpha),Cl(:,i),'k','LineWidth',line_width);
    elseif rem(100*H_array(i),inter)<1e-12
        plot(rad2deg(alpha),Cl(:,i),symbols{floor(i/inter)+1},'LineWidth',line_width);
    end
    hold on;
end

xlabel('$\alpha^{*}, ^{\circ}$','interpreter','latex','FontSize',font_labels);
ylabel('$\overline{C}_{L}$','interpreter','latex','FontSize',font_labels);
xlim([0 max(rad2deg(alpha))]);

Legend = cell(ceil(n_H_array/inter),1);
for iter = 1:n_H_array
    if H_array(iter)==0
        Legend{iter} = 'Steady';
    elseif rem(100*H_array(iter),inter)<1e-12
        Legend{ceil(iter/inter)} = strcat('$H=$', num2str(H_array(iter)));
    end
end

legend(Legend,'Location','best','interpreter','latex','FontSize',font_lgd)
grid on;

% Cd

figure;
colororder(Okabe_Ito)
hold on;
for i = 1:n_H_array
    if i==1
        plot(rad2deg(alpha),Cd(:,i),'k','LineWidth',line_width);
    elseif rem(100*H_array(i),inter)<1e-12
        plot(rad2deg(alpha),Cd(:,i),symbols{floor(i/inter)+1},'LineWidth',line_width);
    end
    hold on;
end

xlabel('$\alpha^{*}, ^{\circ}$','interpreter','latex','FontSize',font_labels);
ylabel('$\overline{C}_{D}$','interpreter','latex','FontSize',font_labels);
xlim([0 max(rad2deg(alpha))]);

legend(Legend,'Location','best','interpreter','latex','FontSize',font_lgd)
grid on;

%% Terms

index = find(H_array==H_display);

% Cl

figure;
colororder(Okabe_Ito)
for i = 1:n_terms
    plot(rad2deg(alpha),term_Cl(:,index,i)-term_Cl(:,1,i),symbols{i},'LineWidth',line_width);
    hold on;
end

xlabel('$\alpha^{*}, ^{\circ}$','interpreter','latex','FontSize',font_labels);
ylabel('$\overline{C}_{L}(\alpha^{*})-C_{L_{s}}(\alpha^{*})$','interpreter','latex','FontSize',font_labels);

legend('$\dot{C}_{v}^{*}\cos\alpha^{*}$','$C_{N_{\alpha}}(\alpha^{*}-Hbk\sin\phi)(1+\sqrt{x_{0}^{*}})^2\cos\alpha^{*}/4$'...
    ,'$\eta C_{N_{\alpha}}(\alpha^{*}-Hbk\sin\phi)^2\sqrt{x_{0}^{*}}\sin\alpha^{*}$'...
    ,'Location','best','interpreter','latex','FontSize',font_lgd)

xlim([0 max(rad2deg(alpha))])
grid on;

% Cd

figure;
colororder(Okabe_Ito)
for i = 1:n_terms
    plot(rad2deg(alpha),term_Cd(:,index,i)-term_Cd(:,1,i),symbols{i},'LineWidth',line_width);
    hold on;
end

xlabel('$\alpha^{*}, ^{\circ}$','interpreter','latex','FontSize',font_labels);
ylabel('$\overline{C}_{D}(\alpha^{*})-C_{D_{s}}(\alpha^{*})$','interpreter','latex','FontSize',font_labels);

legend('$\dot{C}_{v}^{*}\sin\alpha^{*}$','$C_{N_{\alpha}}(\alpha^{*}-Hbk\sin\phi)(1+\sqrt{x_{0}^{*}})^2\sin\alpha^{*}/4$'...
    ,'$-\eta C_{N_{\alpha}}(\alpha^{*}-Hbk\sin\phi)^2\sqrt{x_{0}^{*}}\cos\alpha^{*}$'...
    ,'Location','best','interpreter','latex','FontSize',font_lgd)

xlim([0 max(rad2deg(alpha))])
grid on;