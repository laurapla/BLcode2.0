%% Computing the lift coefficient averaged over a cycle
% University of California, Irvine - Fall 2021
% Laura Pla Olea - lplaolea@uci.edu

clear; clc; close all;

%% Loading files

addpath(genpath('../'),genpath('/data'),genpath('../experimental_data'))
load('dCv.mat');
load('CN_NACA0012_static.mat');

%% Input data

% Geometry
airfoil = ('NACA0012');
M = 0.3; % Mach number
k = 0.1; % Reduced frequency

% Motion
n_A = 11;
A_alpha = linspace(0,10,n_A); % Pitching amplitude [rad]
H = 0; % Plunging amplitude (h/c)
phi = deg2rad(0); % Phase between the pitching and plunging motions [rad]

% Numerical values
[eta, C_Nalpha, ~, ~, ~, alpha1, ~, S1, S2, ~, ~, ~, ~, ~, ~, ~, ~, ~] = input_NACA0012(M);

% Angle parameters
N = 200;
alpha = deg2rad(linspace(0,45,N)).';

%% dCv_matrix processing

n_A_array = length(A_array);
n_k_array = length(k_array);
n_alpha_array = length(alpha_array);

k_index = 1;
for i = 2:n_k_array
    if abs(k-k_array(i))<abs(k-k_array(i-1))
        k_index = i;
    end
end

%% Averaged Cl

% State x10 (separation point)
x10 = separation_point(alpha-H*k*sin(phi),alpha1,S1,S2).';

% State x11 (dCv)
[oldcols,oldrows] = meshgrid(A_array,alpha_array);
[newcols, newrows] = meshgrid(A_alpha,alpha);
x11 = interp2(oldcols,oldrows,dCv_A(:,:,k_index),newcols,newrows);

%% Averaged Cl and Cd

Cl = zeros(N,n_A);
Cd = zeros(N,n_A);

for j = 1:n_A
    
    % Averaged lift coefficient
    Cl(:,j) = x11(:,j).*cos(alpha)...
        +C_Nalpha/2*(alpha-H*k*sin(phi)).*((1+sqrt(x10)).^2/2.*cos(alpha)+2*eta*(alpha-H*k*sin(phi)).*sqrt(x10).*sin(alpha))...
        -deg2rad(A_alpha(j))^2/4*(x11(:,j).*cos(alpha)+8*sin(alpha)/M+C_Nalpha*alpha/2.*(1+sqrt(x10)).^2/2.*cos(alpha)+eta*C_Nalpha*alpha.^2.*sqrt(x10).*sin(alpha))...
        +2*deg2rad(A_alpha(j))*H*k/M*sin(alpha)*sin(phi);

    % Averaged drag coefficient
    Cd(:,j) = x11(:,j).*sin(alpha)...
        +C_Nalpha/2*(alpha-H*k*sin(phi)).*((1+sqrt(x10)).^2/2.*sin(alpha)-eta*(alpha-H*k*sin(phi)).*sqrt(x10).*cos(alpha))...
        -deg2rad(A_alpha(j))^2/4*(x11(:,j).*sin(alpha)-8*cos(alpha)/M+C_Nalpha*alpha/2.*(1+sqrt(x10)).^2/2.*sin(alpha)-eta*C_Nalpha*alpha.^2.*sqrt(x10).*cos(alpha))...
        -2*deg2rad(A_alpha(j))*H*k/M*cos(alpha)*sin(phi);

end

%% Plot Cl

line_width = 1.7;
font_lgd = 10;
font_labels = 14;

symbols = {'-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.'};
Okabe_Ito = [0.902 0.624 0; 0.337 0.737 0.914; 0 0.62 0.451;
    0.941 0.894 0.259; 0 0.447 0.698; 0.835 0.369 0; 0.8 0.475 0.655];

figure;
colororder(Okabe_Ito)
for i = 1:n_A
    plot(rad2deg(alpha),Cl(:,i),symbols{i},'LineWidth',line_width);
    hold on;
end

xlabel('$\alpha^{*}, ^{\circ}$','interpreter','latex','FontSize',font_labels);
ylabel('$\overline{C}_{L}$','interpreter','latex','FontSize',font_labels);

Legend = cell(n_A,1);
for iter = 1:n_A
    Legend{iter} = strcat('$A_{\alpha}=$', num2str((iter-1)), '$^{\circ}$');
end
legend(Legend,'Location','best','interpreter','latex','FontSize',font_lgd)
grid on;

%% Plot Cd

figure;
colororder(Okabe_Ito)
for i = 1:n_A
    hold on;
    plot(rad2deg(alpha),Cd(:,i),symbols{i},'LineWidth',line_width);
end

xlabel('$\alpha^{*}, ^{\circ}$','interpreter','latex','FontSize',font_labels);
ylabel('$\overline{C}_{D}$','interpreter','latex','FontSize',font_labels);

Legend = cell(n_A,1);
for iter = 1:n_A
    Legend{iter} = strcat('$A_{\alpha}=$', num2str((iter-1)), '$^{\circ}$');
end
legend(Legend,'Location','best','interpreter','latex','FontSize',font_lgd)
grid on;