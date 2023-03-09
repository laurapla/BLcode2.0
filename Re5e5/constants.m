%% Data processing to obtain the constants for the Beddoes-Leishman model
% University of California, Irvine - Winter 2023
% Laura Pla Olea - lplaolea@uci.edu

close all; clear; clc;

%% Input

% File with the Cl and Cd data of the airfoil
filename = 'JFM_0012_Re500k.xlsx';

% Convert the file into a Matlab table
[alpha, Cl, Cd] = readvars(filename, 'Range', 'A1:C29');

% Plot Cl and Cd to check they were imported correctly
figure(1);
plot(alpha,Cl,'-o',alpha,Cd,'-o');
xlabel('\alpha, deg'); legend('C_{L}', 'C_{D}', 'Location', 'best');
grid on;

%% Processing

% Normal force coefficient and chord force coefficient
Cn = Cl.*cosd(alpha)+Cd.*sind(alpha);
Cc = Cl.*sind(alpha)-Cd.*cosd(alpha);

% Specify minimum and maximum angles of attack to calculate the slope of
% the linear region
alpha_max = 4; % Maximum angle to calculate the slope [deg]

% Find the index of the angle closest to the specified minimum and maximum
% angles of attack
da = alpha(2)-alpha(1); % Angle increments [deg]
index_alpha_min = find(alpha>=0, 1);
index_alpha_max = find(abs(alpha-alpha_max)<abs(da)/2,1);
alpha_min = alpha(index_alpha_min);
alpha_max = alpha(index_alpha_max);

% Calculate normal force coefficient slope in the linear region
alpha_linear = alpha(index_alpha_min:index_alpha_max);
Cn_linear = Cn(index_alpha_min:index_alpha_max);
Cna_deg = alpha_linear\Cn_linear; % Cn slope [1/deg]
Cna = rad2deg(Cna_deg); % Cn slope [1/rad]
display(Cna);

% Plot Cn and Cc to check if they were calculated correctly
figure(2);
plot(alpha,Cn,'-o',alpha,Cc,'-o',alpha_linear,Cn_linear);
xlabel('\alpha, deg'); legend('C_{N}', 'C_{C}', 'C_{N} linear', 'Location', 'best');
grid on;

%% Separation point

% Separation point according to Kirchhoff's theory
f = (2*sqrt(Cn./(Cna.*deg2rad(alpha)))-1).^2;
f(isnan(f)) = 1;
f(f>1) = 1;
f(f<0) = 0;

% Separation point modeling
alpha1 = interp1(f(3:end),alpha(3:end),0.7);
index_alpha1 = find(abs(alpha-alpha1)<abs(da)/2,1);
if f(index_alpha1)<0.7
    index_alpha1 = index_alpha1-1;
end
ft = fittype(@(a1,S1,S2,x) piecewiseLine(x,a1,S1,S2));
fun = fit(alpha, f, ft, 'StartPoint', [14, 3, 1] );
alpha1 = fun.a1;
S1 = fun.S1;
S2 = fun.S2;
display(alpha1);
display(S1);
display(S2);

% Reconstruct the separation point modeling
n = length(alpha);
f_model = zeros(n,1);
for i = 1:n
    if alpha(i)<alpha1
        f_model(i) = 1-0.3*exp((abs(alpha(i))-alpha1)/S1);
    else
        f_model(i) = 0.04+0.66*exp(alpha1-(abs(alpha(i)))/S2);
    end
end

% Plot separation point to check if it was calculated and modeled correctly
figure(3);
plot(alpha,f,'-o',alpha,f_model);
xlabel('\alpha, deg'); legend('f', 'f model', 'Location', 'best');
grid on;