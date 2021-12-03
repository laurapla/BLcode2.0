clear; clc; close all;

% Data Processing
interval = 5;
number = 50;

figure;
for i = 0:number
    
    filename = ['data/dCv_k' num2str(i)];
    load(filename);
    
    N = length(alpha_array);
    n_t = length(t);
    
    if rem(i,interval)==0
        hold on;
        if i==0
            plot(rad2deg(alpha_array),zeros(1,N))
        else
            plot(rad2deg(alpha_array),rms(fft(dCv_array)))
        end
    end
    
end