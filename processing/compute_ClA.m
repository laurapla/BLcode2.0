%% Computing the lift coefficient averaged over a cycle
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

Cl = zeros(N,n_A);

for j = 1:n_A
    
    % Averaged lift coefficient
    Cl(:,j) = x11(:,j).*cos(alpha)...
            +C_Nalpha/2*(alpha-H*k*sin(phi)).*((1+sqrt(x10)).^2/2.*cos(alpha)+2*eta*(alpha-H*k*sin(phi)).*sqrt(x10).*sin(alpha))...
            -deg2rad(A_alpha(j))^2/4*(x11(:,j).*cos(alpha)+sin(alpha)/(8*M)+C_Nalpha*alpha/2.*(1+sqrt(x10)).^2/2.*cos(alpha)+eta*C_Nalpha*alpha.^2.*sqrt(x10).*sin(alpha))...
            +deg2rad(A_alpha(j))*H*k/M*sin(alpha)*sin(phi);
    
end

%% Static value

alpha_s = CN_NACA0012_static(:,1);
Cl_s = CN_NACA0012_static(:,2);

figure;
plot(alpha_s,Cl_s,'--o');
 for i = 1:n_A
        hold on;
        plot(rad2deg(alpha),Cl(:,i));
 end
 
 xlabel('\alpha^{*}'); ylabel('$\overline{C}_{N}$','interpreter','latex');
 
 Legend = cell(n_A+1,1);
for iter = 1:n_A+1
    if iter==1
        Legend{iter} = strcat('Steady Cn');
    else
        Legend{iter} = strcat('A_{\alpha}=', num2str((iter-2)), 'º');
    end
 end
 legend(Legend,'Location','bestoutside')