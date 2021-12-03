% Finding the mean value of dCv over a cycle
% University of California, Irvine - Fall 2021
% Laura Pla Olea - lplaolea@uci.edu

clear; clc; close all;

% Parameters of Data Processing
parameter = 'H';
if strcmp(parameter,'A')
    interval = 2;
    number = 15;
elseif strcmp(parameter,'k')
    interval = 5;
    number = 50;
elseif strcmp(parameter,'H')
    interval = 1;
    number = 11;
end

% Storing the results
dCv = zeros(200,number);
para = zeros(1,number);

figure(1);
for i = 0:number
    
    if strcmp(parameter,'A')
        filename = ['data/dCv_A' num2str(i)];
    elseif strcmp(parameter,'k')
        filename = ['data/dCv_k' num2str(i)];
    elseif strcmp(parameter,'H')
        filename = ['data/dCv_H' num2str(5*i)];
    end
    load(filename);
    
    N = length(alpha_array);
    n_t = length(t);
    
    % Storing the results
    dCv(:,i+1) = mean(dCv_array);
    dCv(isinf(dCv)|isnan(dCv)) = 0;
    if strcmp(parameter,'A')
        para(i+1) = i;
    elseif strcmp(parameter,'k')
        para(i+1) = i/100;
    elseif strcmp(parameter,'H')
        para(i+1) = i*5/100;
    end
    
    if rem(i,interval)==0
        hold on;
        plot(rad2deg(alpha_array),dCv(:,i+1))
    end
    
end

%% Plot parameters

Legend = cell(floor(number/interval)+1,1);
for iter=1:number+1
    if rem(iter-1,interval)==0
        if strcmp(parameter,'A')
            Legend{(iter-1)/interval+1} = strcat('A_{\alpha}=', num2str((iter-1)), 'º');
        elseif strcmp(parameter,'k')
            Legend{(iter-1)/interval+1} = strcat('k=', num2str((iter-1)/100));
        elseif strcmp(parameter,'H')
            Legend{(iter-1)/interval+1} = strcat('H=', num2str(5*(iter-1)/100));
        end
    end
 end
 legend(Legend,'Location','best')
 
 xlabel('\alpha [º]'); ylabel('$\dot{C}_{v}^{*}$','interpreter','latex');
 
 %% Plot as a function of the parameter
 
 interval_angle = 15;
 N_angle = N;
 fffit = zeros(3,floor(N_angle/interval_angle));
 
 figure(2);
 for i = 1:N_angle
     f = fit(para.',dCv(i,:).','poly2');
     fffit(:,i) = coeffvalues(f);
     if rem(i-1,interval_angle)==0
        hold on;
        plot(para,dCv(i,:));
     end
 end
 
 Legend2 = cell(floor(N_angle/interval_angle),1);
 for iter = 1:N_angle
     if rem(iter-1,interval_angle)==0
        Legend2{(iter-1)/interval_angle+1} = strcat('\alpha^{*}=', num2str(round(rad2deg(alpha_array(iter)))), 'º');
     end
 end
 legend(Legend2,'Location','bestoutside')
 
 if strcmp(parameter,'A')
     xlabel('A_{\alpha} [º]');
 elseif strcmp(parameter,'k')
     xlabel('k');
 elseif strcmp(parameter,'H')
     xlabel('H');
     xlim([0 .55])
 end
 
 ylabel('$\dot{C}_{v}^{*}$','interpreter','latex');