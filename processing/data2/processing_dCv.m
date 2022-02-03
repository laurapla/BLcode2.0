% Processing the results of dCv as a function of k, A_alpha, and H
% University of California, Irvine - Winter 2022
% Laura Pla Olea - lplaolea@uci.edu

clear; clc; close all
addpath(genpath('../'),genpath('/data2'))

%% Get and store data

n_k = 11;
k_array = linspace(0,0.5,n_k);

n_alpha = 20;
alpha_array = linspace(0,pi/4,n_alpha);

n_A = 11;
dCv_A_array = zeros(n_alpha,n_A);
dCv_A = zeros(n_alpha,n_A,n_k);

n_H = 12;
dCv_H_array = zeros(n_alpha,n_H);
dCv_H = zeros(n_alpha,n_H,n_k);

for i = 1:n_k
    
    filenameA = ['dCv_A_k',num2str(5*(i-1))];
    load(filenameA);
    
    % Storing the results
    for j = 1:n_A
        
        dCv_array = dCv_matrix(:,:,j);
        
        % Storing the results
        dCv_A_array(:,j) = mean(dCv_array);
        dCv_A_array(isinf(dCv_A_array)|isnan(dCv_A_array)) = 0;
        
    end
    
    dCv_A(:,:,i) = dCv_A_array;
    
    
    filenameH = ['dCv_H_k',num2str(5*(i-1))];
    load(filenameH);
    
    % Storing the results
    for j = 1:n_H
        
        dCv_array = dCv_matrix(:,:,j);
        
        % Storing the results
        dCv_H_array(:,j) = mean(dCv_array);
        dCv_H_array(isinf(dCv_H_array)|isnan(dCv_H_array)) = 0;
        
    end
    
    dCv_H(:,:,i) = dCv_H_array;
    
end

save('dCv.mat','alpha_array','k_array','A_array','H_array','dCv_A','dCv_H');

%% Plot dCv vs. alpha as a function of k (for a given A_alpha)

A_index = 3; % Determines the A_alpha for which we plot the results

Legend_alpha = cell(n_alpha,1);

figure(1);
for i = 1:n_alpha
    
    plot(k_array,squeeze(dCv_A(i,A_index,:)));
    hold on;
    Legend_alpha{i} = strcat('\alpha^{*}=', num2str(round(rad2deg(alpha_array(i)))), 'º');
    
end

legend(Legend_alpha,'Location','bestoutside')
title(['A_{\alpha}=',num2str(A_array(A_index)),'º, H=0']);
xlabel('k'); ylabel('$\dot{C}_{v}^{*}$','interpreter','latex');

%% Plot dCv vs. alpha as a function of k (for a given H)

H_index = 3; % Determines the H for which we plot the results

figure(2);
for i = 1:n_alpha
    
    plot(k_array,squeeze(dCv_H(i,H_index,:)));
    hold on;
    
end

legend(Legend_alpha,'Location','bestoutside')
title(['A_{\alpha}=0º, H=',num2str(H_array(H_index))]);
xlabel('k'); ylabel('$\dot{C}_{v}^{*}$','interpreter','latex');

%% Plot dCv vs. alpha as a function of A_alpha (for a given k)

k_index = 3; % Determines the k for which we plot the results

Legend_A = cell(n_A,1);

figure(3);
for i = 1:n_A
    
    plot(alpha_array,squeeze(dCv_A(:,i,k_index)));
    hold on;
    Legend_A{i} = strcat('A_{\alpha}=', num2str(round(A_array(i))), 'º');
    
end

legend(Legend_A,'Location','best')
title(['k=',num2str(k_array(k_index)),', H=0']);
xlabel('\alpha^{*}'); ylabel('$\dot{C}_{v}^{*}$','interpreter','latex');

%% Plot dCv as a function of the plunging amplitude (for a given k)

k_index = 7; % Determines the k for which we plot the results

Legend_H = cell(n_H,1);

figure(4);
for i = 1:n_H
    
    plot(alpha_array,squeeze(dCv_H(:,i,k_index)));
    hold on;
    Legend_H{i} = strcat('H=', num2str(H_array(i)));
    
end

legend(Legend_H,'Location','best')
title(['k=',num2str(k_array(k_index)),', A_{\alpha}=0']);
xlabel('\alpha^{*}'); ylabel('$\dot{C}_{v}^{*}$','interpreter','latex');

%% Plot dCv vs. A_alpha as a function of alpha (for a given k)

k_index = 3; % Determines the k for which we plot the results

figure(5);
for i = 1:n_alpha
    
    plot(A_array,squeeze(dCv_A(i,:,k_index)));
    hold on;
    
end

legend(Legend_alpha,'Location','bestoutside')
title(['k=',num2str(k_array(k_index)),', H=0']);
xlabel('A_{\alpha}'); ylabel('$\dot{C}_{v}^{*}$','interpreter','latex');

%% Plot dCv vs. H as a function of alpha (for a given k)

k_index = 7; % Determines the k for which we plot the results

figure(6);
for i = 1:n_alpha
    
    plot(H_array,squeeze(dCv_H(i,:,k_index)));
    hold on;
    
end

legend(Legend_alpha,'Location','bestoutside')
title(['k=',num2str(k_array(k_index)),', A_{\alpha}=0º']);
xlabel('H'); ylabel('$\dot{C}_{v}^{*}$','interpreter','latex');