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
alpha = deg2rad(linspace(0,45,N));

%% Averaged Cl

% State x10 (separation point)
x10 = separation_point(alpha-H*k*sin(phi),alpha1,S1,S2);

% State x11 (dCv)
[oldcols,oldrows] = meshgrid(A_array,alpha_array);
[newcols, newrows] = meshgrid(A_alpha,alpha);
x11 = interp2(oldcols,oldrows,dCv_alphaA,newcols,newrows);

%% Averaged Cl

Cl = zeros(N,n_A);

for j = 1:n_A
    for i = 1:N
        
        % Averaged lift coefficient
        Cl(i,j) = x11(i,j)*cos(alpha(i))...
            +C_Nalpha/2*(alpha(i)-H*k*sin(phi))*((1+sqrt(x10(i)))^2/2*cos(alpha(i))+2*eta*(alpha(i)-H*k*sin(phi))*sqrt(x10(i))*sin(alpha(i)))...
            -deg2rad(A_alpha(j))^2/4*(x11(i,j)*cos(alpha(i))+sin(alpha(i))/(8*M)+C_Nalpha*alpha(i)/2*(1+sqrt(x10(i)))^2/2*cos(alpha(i))+eta*C_Nalpha*alpha(i)^2*sqrt(x10(i))*sin(alpha(i)))...
            +deg2rad(A_alpha(j))*H*k/M*sin(alpha(i))*sin(phi);
        
    end
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
 
 xlabel('\alpha^{*}'); ylabel('$\overline{C}_{L}$','interpreter','latex');
 
 Legend = cell(n_A+1,1);
for iter = 1:n_A+1
    if iter==1
        Legend{iter} = strcat('Steady Cl');
    else
        Legend{iter} = strcat('A_{\alpha}=', num2str((iter-2)), 'º');
    end
 end
 legend(Legend,'Location','bestoutside')