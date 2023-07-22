%% Computing the lift coefficient averaged over a cycle
% University of California, Irvine - Fall 2021
% Laura Pla Olea - lplaolea@uci.edu

clear; clc; close all;

%% Loading files

addpath(genpath('../'),genpath('/data'),genpath('../experimental_data'))
load('dCv_ak.mat');
load('CN_NACA0012_static.mat');

%% Input data

% Geometry
airfoil = ('NACA0012');
M = 0.3; % Mach number
n_k = 51;
k = linspace(0,0.5,n_k); % Reduced frequency

% Motion
A_alpha = 5; % Pitching amplitude [rad]
H = 0; % Plunging amplitude (h/c)
phi = deg2rad(0); % Phase between the pitching and plunging motions [rad]

% Numerical values
[eta, C_Nalpha, ~, ~, ~, alpha1, ~, S1, S2, ~, ~, ~, ~, ~, ~, ~, ~, ~] = input_NACA0012(M);

% Angle parameters
N = 200;
alpha = deg2rad(linspace(0,45,N));

%% Averaged Cl

% State x11 (dCv)
[oldcols,oldrows] = meshgrid(k_array,alpha_array);
[newcols, newrows] = meshgrid(k,alpha);
x11 = interp2(oldcols,oldrows,dCv_alphak,newcols,newrows);

%% Averaged Cl

Cl = zeros(N,n_k);
x10 = zeros(1,N);

for j = 1:n_k
    for i = 1:N
        
        % State x10 (separation point)
        x10(i) = separation_point(alpha(i)-H*k(j)*sin(phi),alpha1,S1,S2);
        
        % Averaged lift coefficient
        Cl(i,j) = x11(i,j)*cos(alpha(i))...
            +C_Nalpha/2*(alpha(i)-H*k(j)*sin(phi))*((1+sqrt(x10(i)))^2/2*cos(alpha(i))+2*eta*(alpha(i)-H*k(j)*sin(phi))*sqrt(x10(i))*sin(alpha(i)))...
            -deg2rad(A_alpha)^2/4*(x11(i,j)*cos(alpha(i))+sin(alpha(i))/(8*M)+C_Nalpha*alpha(i)/2*(1+sqrt(x10(i)))^2/2*cos(alpha(i))+eta*C_Nalpha*alpha(i)^2*sqrt(x10(i))*sin(alpha(i)))...
            +deg2rad(A_alpha)*H*k(j)/M*sin(alpha(i))*sin(phi);
        
    end
end

%% Static value

alpha_s = CN_NACA0012_static(:,1);
Cl_s = CN_NACA0012_static(:,2);

figure;
plot(alpha_s,Cl_s);
 for i = 1:n_k
        hold on;
        plot(rad2deg(alpha),Cl(:,i));
 end
 
 xlabel('\alpha^{*}'); ylabel('$\overline{C}_{L}$','interpreter','latex');
 
 Legend = cell(n_k+1,1);
for iter = 1:n_k+1
    if iter==1
        Legend{iter} = strcat('Steady Cl');
    else
        Legend{iter} = strcat('k=', num2str((iter-2)/100));
    end
 end
 legend(Legend,'Location','bestoutside')