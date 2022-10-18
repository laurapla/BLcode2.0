%% Computing the lift coefficient averaged over a cycle
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

n_H_array = length(H_array);

% Motion
A_alpha = 0; % Pitching amplitude [rad]
H = linspace(0,0.55,n_H_array); % Plunging amplitude (h/c)
phi = deg2rad(0); % Phase between the pitching and plunging motions [rad]

% Numerical values
[eta, C_Nalpha, ~, ~, ~, alpha1, ~, S1, S2, ~, ~, ~, ~, ~, ~, ~, ~, ~] = input_NACA0012(M);

% Angle parameters
N = 200;
alpha = deg2rad(linspace(0,45,N));

%% dCv_matrix processing

n_k_array = length(k_array);

k_index = 1;
for i = 2:n_k_array
    if abs(k-k_array(i))<abs(k-k_array(i-1))
        k_index = i;
    end
end

%% State x11 (dCv)

[oldcols,oldrows] = meshgrid(H_array,alpha_array);
[newcols, newrows] = meshgrid(H,alpha);
x11 = interp2(oldcols,oldrows,dCv_H(:,:,k_index),newcols,newrows);

%% Average Cl

x10 = zeros(1,N);
Cl = zeros(N,n_H_array);
Cd = zeros(N,n_H_array);

for j = 1:n_H_array
    for i = 1:N
        
        % State x10 (separation point)
        x10(i) = separation_point(alpha(i)-H(j)*k*sin(phi),alpha1,S1,S2);
        
        % Averaged lift coefficient
        Cl(i,j) = x11(i,j)*cos(alpha(i))...
            +C_Nalpha/2*(alpha(i)-H(j)*k*sin(phi))*((1+sqrt(x10(i)))^2/2*cos(alpha(i))+2*eta*(alpha(i)-H(j)*k*sin(phi))*sqrt(x10(i))*sin(alpha(i)))...
            -deg2rad(A_alpha)^2/4*(x11(i,j)*cos(alpha(i))+sin(alpha(i))/(8*M)+C_Nalpha*alpha(i)/2*(1+sqrt(x10(i)))^2/2*cos(alpha(i))+eta*C_Nalpha*alpha(i)^2*sqrt(x10(i))*sin(alpha(i)))...
            +deg2rad(A_alpha)*H(j)*k/M*sin(alpha(i))*sin(phi);

        % Averaged drag coefficient
        Cd(i,j) = x11(i,j).*sin(alpha(i))...
            +C_Nalpha/2*(alpha(i)-H(j)*k*sin(phi)).*((1+sqrt(x10(i))).^2/2.*sin(alpha(i))-eta*(alpha(i)-H(j)*k*sin(phi)).*sqrt(x10(i)).*cos(alpha(i)))...
            -deg2rad(A_alpha)^2/4*(x11(i,j).*sin(alpha(i))-8*cos(alpha(i))/M+C_Nalpha*alpha(i)/2.*(1+sqrt(x10(i))).^2/2.*sin(alpha(i))-eta*C_Nalpha*alpha(i).^2.*sqrt(x10(i)).*cos(alpha(i)))...
            -2*deg2rad(A_alpha)*H(j)*k/M*cos(alpha(i))*sin(phi);

    end
end

%% Plot Cl

line_width = 1.7;
font_lgd = 10;
font_labels = 14;
inter = 2;

symbols = {'-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.'};
Okabe_Ito = [0.902 0.624 0; 0.337 0.737 0.914; 0 0.62 0.451;
    0.941 0.894 0.259; 0 0.447 0.698; 0.835 0.369 0; 0.8 0.475 0.655];

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

xlabel('$\alpha^{*}$','interpreter','latex','FontSize',font_labels);
ylabel('$\overline{C}_{L}$','interpreter','latex','FontSize',font_labels);

Legend = cell(ceil(n_H_array/inter),1);
for iter = 1:n_H_array
    if H_array(iter)==0
        Legend{iter} = 'Steady';
    elseif rem(100*H_array(iter),inter)<1e-12
        Legend{ceil(iter/inter)} = strcat('$H=$', num2str(5*(iter-1)/100));
    end
end
legend(Legend,'Location','best','interpreter','latex','FontSize',font_lgd)
grid on;

%% Plot Cd

figure;
colororder(Okabe_Ito)
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

legend(Legend,'Location','best','interpreter','latex','FontSize',font_lgd)
grid on;