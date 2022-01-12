%% Simulation of the Beddoes-Leishman State Space Model (NACA0012)
% University of California, Irvine - Fall 2021
% Laura Pla Olea - lplaolea@uci.edu

clear; clc; close all;

%% Data

load('data/dCv_A10');

N = length(alpha_array);
n_t = length(t)/10;
Fs = 1/(t(2)-t(1));

%% FFT

Y = fft(dCv_array(:,1));

P2 = abs(Y/N);
P1 = P2(1:n_t/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(n_t/2))/n_t;
figure; plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% FFT 2

X = 1/N*fftshift(fft(dCv_array(:,1),N));%N-point complex DFT

df=Fs/N; %frequency resolution
sampleIndex = -N/2:N/2-1; %ordered index for FFT plot
foo=sampleIndex*df; %x-axis index converted to ordered frequencies
figure; stem(foo,abs(X)); %magnitudes vs frequencies
xlabel('f (Hz)'); ylabel('|X(k)|'); xlim([-foo(end) foo(end)])