%% Computing the lift coefficient averaged over a cycle
% University of California, Irvine - Fall 2021
% Laura Pla Olea - lplaolea@uci.edu

clear; clc; close all;

%% Loading files

addpath(genpath('../'),genpath('/data'),genpath('../experimental_data'))
load('dCv_H.mat');
load('CN_NACA0012_static.mat');

%% Input data

% Geometry
airfoil = ('NACA0012');
M = 0.3; % Mach number
k = 0.5; % Reduced frequency

% Motion
n_H = 12;
A_alpha = 0; % Pitching amplitude [rad]
H = linspace(0,0.55,n_H); % Plunging amplitude (h/c)
phi = deg2rad(0); % Phase between the pitching and plunging motions [rad]

% Numerical values
[eta, C_Nalpha, ~, ~, ~, alpha1, ~, S1, S2, ~, ~, ~, ~, ~, ~, ~, ~, ~] = input_NACA0012(M);

% Angle parameters
N = 200;
alpha = deg2rad(linspace(0,45,N));

%% Averaged Cl

% State x11 (dCv)
[oldcols,oldrows] = meshgrid(H_array,alpha_array);
[newcols, newrows] = meshgrid(H,alpha);
x11 = interp2(oldcols,oldrows,dCv,newcols,newrows);

%% Averaged Cl

Cl = zeros(N,n_H);
x10 = zeros(1,N);

for j = 1:n_H
    for i = 1:N
        
        % State x10 (separation point)
        x10(i) = separation_point(alpha(i)-H(j)*k*sin(phi),alpha1,S1,S2);
        
        % Averaged lift coefficient
        Cl(i,j) = x11(i,j)*cos(alpha(i))...
            +C_Nalpha/2*(alpha(i)-H(j)*k*sin(phi))*((1+sqrt(x10(i)))^2/2*cos(alpha(i))+2*eta*(alpha(i)-H(j)*k*sin(phi))*sqrt(x10(i))*sin(alpha(i)))...
            -deg2rad(A_alpha)^2/4*(x11(i,j)*cos(alpha(i))+sin(alpha(i))/(8*M)+C_Nalpha*alpha(i)/2*(1+sqrt(x10(i)))^2/2*cos(alpha(i))+eta*C_Nalpha*alpha(i)^2*sqrt(x10(i))*sin(alpha(i)))...
            +deg2rad(A_alpha)*H(j)*k/M*sin(alpha(i))*sin(phi);
        
    end
end

%% Static value

alpha_s = CN_NACA0012_static(:,1);
Cl_s = CN_NACA0012_static(:,2);

figure;
plot(alpha_s,Cl_s);
for i = 1:n_H
    hold on;
    plot(rad2deg(alpha),Cl(:,i));
end

xlabel('\alpha^{*}'); ylabel('$\overline{C}_{L}$','interpreter','latex');

Legend = cell(n_H+1,1);
for iter = 1:n_H+1
    if iter==1
        Legend{iter} = strcat('Steady Cl');
    else
        Legend{iter} = strcat('H=', num2str(5*(iter-2)/100));
    end
end
legend(Legend,'Location','bestoutside')
xlim([0 25])
ylim([0 1.6])

%%

Cl35 = [Cl(1:74,8); Cl(76:end,8)]; a35 = [alpha(1:74) alpha(76:end)];
Cl50 = [Cl(1:75,11); Cl(79:end,11)]; a50 = [alpha(1:75) alpha(79:end)];

figure;
plot(alpha_s,Cl_s,'--o');
hold on; plot(rad2deg(alpha),Cl(:,3));
hold on; plot(rad2deg(a35),Cl35);
hold on; plot(rad2deg(a50),Cl50);
xlim([0 25])
xlabel('\alpha^*'); ylabel('$\overline{C}_{L}$','interpreter','latex');
grid on
legend('Steady C_{L}','H=0.1','H=0.35','H=0.5','Location','southeast')