% Processing the results of dCv as a function of k, A_alpha, and H
% University of California, Irvine - Winter 2022
% Laura Pla Olea - lplaolea@uci.edu

clear; clc; close all
addpath(genpath('../'),genpath('/data/M0.7'))

%% Get and store data

n_k = 21;
k_array = linspace(0,1,n_k);

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

width = 1.7;
font_lgd = 10;
font_labels = 14;

symbols = {'-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.','-', '--', ':','-.'};
Okabe_Ito = [0.902 0.624 0; 0.337 0.737 0.914; 0 0.62 0.451;
0.941 0.894 0.259; 0 0.447 0.698; 0.835 0.369 0; 0.8 0.475 0.655];

A_index = 3; % Determines the A_alpha for which we plot the results

Legend_alpha = cell(12,1);

figure(1);
colororder(Okabe_Ito)
for i = 1:12

    if i==1
        plot(k_array,squeeze(dCv_A(i,A_index,:)),'k','LineWidth',width);
    else
        plot(k_array,squeeze(dCv_A(i,A_index,:)),symbols{i},'LineWidth',width);
    end
    hold on;
    Legend_alpha{i} = strcat('$\alpha^{*}=$', num2str(round(rad2deg(alpha_array(i)))), '$^{\circ}$');
    
end

legend(Legend_alpha,'Location','best','interpreter','latex','FontSize',font_lgd)
title(['A_{\alpha}=',num2str(A_array(A_index)),'º, H=0']);
xlabel('$k$','interpreter','latex','FontSize',font_labels);
ylabel('$\dot{C}_{v}^{*}$','interpreter','latex','FontSize',font_labels);
grid on;

%% Plot dCv vs. alpha as a function of k (for a given H)

H_index = 3; % Determines the H for which we plot the results

figure(2);
colororder(Okabe_Ito)
for i = 1:12
    
    if i==1
        plot(k_array,squeeze(dCv_H(i,H_index,:)),'k','LineWidth',width);
    else
        plot(k_array,squeeze(dCv_H(i,H_index,:)),symbols{i},'LineWidth',width);
    end
    hold on;
    
end

legend(Legend_alpha,'Location','best','interpreter','latex','FontSize',font_lgd)
title(['A_{\alpha}=0º, H=',num2str(H_array(H_index))]);
xlabel('$k$','interpreter','latex','FontSize',font_labels);
ylabel('$\dot{C}_{v}^{*}$','interpreter','latex','FontSize',font_labels);
grid on;

%% Plot dCv vs. alpha as a function of A_alpha (for a given k)

k_index = 3; % Determines the k for which we plot the results

Legend_A = cell(7,1);

figure(3);
colororder(Okabe_Ito)
for i = 1:7
    
    if i==1
        plot(rad2deg(alpha_array),squeeze(dCv_A(:,i,k_index)),'k','LineWidth',width);
    else
        plot(rad2deg(alpha_array),squeeze(dCv_A(:,i,k_index)),symbols{i},'LineWidth',width);
    end
    hold on;
    Legend_A{i} = strcat('$A_{\alpha}=$', num2str(round(A_array(i))), '$^{\circ}$');
    
end

legend(Legend_A,'Location','best','interpreter','latex','FontSize',font_lgd)
title(['k=',num2str(k_array(k_index)),', H=0']);
xlabel('$\alpha^{*}, ^{\circ}$','interpreter','latex','FontSize',font_labels);
ylabel('$\dot{C}_{v}^{*}$','interpreter','latex','FontSize',font_labels);
grid on;

%% Plot dCv as a function of the plunging amplitude (for a given k)

k_index = 7; % Determines the k for which we plot the results

Legend_H = cell(n_H,1);

figure(4);
colororder(Okabe_Ito)
for i = 1:n_H
    
    if i==1
        plot(rad2deg(alpha_array),squeeze(dCv_H(:,i,k_index)),'k','LineWidth',width);
    else
        plot(rad2deg(alpha_array),squeeze(dCv_H(:,i,k_index)),symbols{i},'LineWidth',width);
    end
    hold on;
    Legend_H{i} = strcat('$H=$', num2str(H_array(i)));
    
end

legend(Legend_H,'Location','best','interpreter','latex','FontSize',font_lgd)
title(['k=',num2str(k_array(k_index)),', A_{\alpha}=0']);
xlabel('$\alpha^{*}, ^{\circ}$','interpreter','latex','FontSize',font_labels);
ylabel('$\dot{C}_{v}^{*}$','interpreter','latex','FontSize',font_labels);
grid on;

%% Plot dCv vs. A_alpha as a function of alpha (for a given k)

k_index = 3; % Determines the k for which we plot the results

figure(5);
colororder(Okabe_Ito)
for i = 1:12
    
    if i==1
        plot(A_array,squeeze(dCv_A(i,:,k_index)),'k','LineWidth',width);
    else
        plot(A_array,squeeze(dCv_A(i,:,k_index)),symbols{i},'LineWidth',width);
    end
    hold on;
    
end

legend(Legend_alpha,'Location','best','interpreter','latex','FontSize',font_lgd)
title(['k=',num2str(k_array(k_index)),', H=0']);
xlabel('$A_{\alpha}$','interpreter','latex','FontSize',font_labels);
ylabel('$\dot{C}_{v}^{*}$','interpreter','latex','FontSize',font_labels);
grid on;

%% Plot dCv vs. H as a function of alpha (for a given k)

k_index = 7; % Determines the k for which we plot the results

figure(6);
colororder(Okabe_Ito)
for i = 1:12
    
    if i==1
        plot(H_array,squeeze(dCv_H(i,:,k_index)),'k','LineWidth',width);
    else
        plot(H_array,squeeze(dCv_H(i,:,k_index)),symbols{i},'LineWidth',width);
    end
    hold on;
    
end

legend(Legend_alpha,'Location','best','interpreter','latex','FontSize',font_lgd)
title(['k=',num2str(k_array(k_index)),', A_{\alpha}=0º']);
xlabel('$H$','interpreter','latex','FontSize',font_labels);
ylabel('$\dot{C}_{v}^{*}$','interpreter','latex','FontSize',font_labels);
xlim([H_array(1) H_array(end)])
grid on;