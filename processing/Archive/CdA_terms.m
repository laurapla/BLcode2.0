%% Separately computing the terms of the drag coefficient averaged over a cycle
% University of California, Irvine - Fall 2022
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

%% Pre-calculations

n_A_array = length(A_array);
n_k_array = length(k_array);
n_alpha_array = length(alpha_array);

% Motion
A_alpha = linspace(0,10,n_A_array); % Pitching amplitude [rad]
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
x10 = separation_point(alpha,alpha1,S1,S2).';

% State x11 (dCv)
[oldcols,oldrows] = meshgrid(A_array,alpha_array);
[newcols, newrows] = meshgrid(A_alpha,alpha);
x11 = interp2(oldcols,oldrows,dCv_A(:,:,k_index),newcols,newrows);

%% Average Cd

n_terms = 7;
term_Cd = zeros(N,n_A_array,n_terms);

for j = 1:n_A_array
    
    % Average drag coefficient
    term_Cd(:,j,1) = x11(:,j).*sin(alpha);
    term_Cd(:,j,2) = C_Nalpha/2*alpha.*(1+sqrt(x10)).^2/2.*sin(alpha);
    term_Cd(:,j,3) = -C_Nalpha/2*alpha*2*eta.*alpha.*sqrt(x10).*cos(alpha);
    term_Cd(:,j,4) = -deg2rad(A_alpha(j))^2/4*x11(:,j).*sin(alpha);
    term_Cd(:,j,5) = deg2rad(A_alpha(j))^2/4*8*cos(alpha)/M;
    term_Cd(:,j,6) = -deg2rad(A_alpha(j))^2/4*C_Nalpha*alpha/2.*(1+sqrt(x10)).^2/2.*sin(alpha);
    term_Cd(:,j,7) = deg2rad(A_alpha(j))^2/4*eta*C_Nalpha*alpha.^2.*sqrt(x10).*cos(alpha);
    
end

Cd = zeros(N,n_A_array);
for i = 1:n_terms
    Cd = Cd+term_Cd(:,:,i);
end

%% Static value

line_width = 1.7;
font_lgd = 10;
font_labels = 14;
inter = 2;

symbols = {'-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.'};
Okabe_Ito = [0.902 0.624 0; 0.337 0.737 0.914; 0 0.62 0.451;
    0.941 0.894 0.259; 0 0.447 0.698; 0.835 0.369 0; 0.8 0.475 0.655];

figure;
colororder(Okabe_Ito)
hold on;
for i = 1:n_A_array
    if i==1
        plot(rad2deg(alpha),Cd(:,i),'k','LineWidth',line_width);
    elseif rem(A_array(i),inter)==0
        plot(rad2deg(alpha),Cd(:,i),symbols{floor(i/inter)+1},'LineWidth',line_width);
    end
    hold on;
end

xlabel('$\alpha^{*}, ^{\circ}$','interpreter','latex','FontSize',font_labels);
ylabel('$\overline{C}_{D}$','interpreter','latex','FontSize',font_labels);
xlim([0 max(rad2deg(alpha))]);

Legend = cell(ceil(n_A_array/inter),1);
for iter = 1:n_A_array
    if A_array(iter)==0
        Legend{iter} = 'Steady';
    elseif rem(A_array(iter),inter)==0
        Legend{ceil(iter/inter)} = strcat('$A_{\alpha}=$', num2str(iter-1), '$^{\circ}$');
    end
end
legend(Legend,'Location','best','interpreter','latex','FontSize',font_lgd)
grid on;

%% Terms

index = 6;

figure;
colororder(Okabe_Ito)
for i = 1:n_terms
    plot(rad2deg(alpha),term_Cd(:,index,i),symbols{i},'LineWidth',line_width);
    hold on;
end

xlabel('$\alpha^{*}, ^{\circ}$','interpreter','latex','FontSize',font_labels);

legend('$\dot{C}_{v}^{*}\sin\alpha^{*}$','$C_{N_{\alpha}}\alpha^{*}(1+\sqrt{x_{0}^{*}})^2\sin\alpha^{*}/4$'...
    ,'$-\eta C_{N_{\alpha}}\alpha^{*2}\sqrt{x_{0}^{*}}\cos\alpha^{*}$'...
    ,'$-A_{\alpha}^{2}\dot{C}_{v}^{*}\sin\alpha^{*}/4$','$2A_{\alpha}^{2}\cos\alpha^{*}/M$'...
    ,'$-A_{\alpha}^{2}C_{N_{\alpha}}\alpha^{*}(1+\sqrt{x_{0}^{*}})^2\sin\alpha^{*}/16$'...
    ,'$A_{\alpha}^{2}\eta C_{N_{\alpha}}\alpha^{*2}\sqrt{x_{0}^{*}}\cos\alpha^{*}/4$'...
    ,'Location','best','interpreter','latex','FontSize',font_lgd,'NumColumns',2)

xlim([0 max(rad2deg(alpha))])
grid on;

%% Grouping the A_alpha terms

figure;
colororder(Okabe_Ito)
for i = 1:4
    if i==4
        plot(rad2deg(alpha),term_Cd(:,index,4)+term_Cd(:,index,5)+term_Cd(:,index,6)+term_Cd(:,index,7),symbols{i},'LineWidth',line_width);
    else
        plot(rad2deg(alpha),term_Cd(:,index,i),symbols{i},'LineWidth',line_width);
    end
    hold on;
end

xlabel('$\alpha^{*}, ^{\circ}$','interpreter','latex','FontSize',font_labels);

legend('$\dot{C}_{v}^{*}\sin\alpha^{*}$','$C_{N_{\alpha}}\alpha^{*}(1+\sqrt{x_{0}^{*}})^2\sin\alpha^{*}/4$'...
    ,'$-\eta C_{N_{\alpha}}\alpha^{*2}\sqrt{x_{0}^{*}}\cos\alpha^{*}$'...
    ,'$A_{\alpha}^{2}/4\times\chi_{DA_{\alpha}}$'...
    ,'Location','best','interpreter','latex','FontSize',font_lgd)

xlim([0 max(rad2deg(alpha))])
grid on;