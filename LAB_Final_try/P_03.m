clear all; clc; close all;
load('data_given');

signal_freq = 1e6;
DMA_freq = 16e6;
DAC_freq = 128e6;
Fc = 32e6;
Fif = 2e6;

M_1 = DMA_freq/signal_freq;
M_2 = DAC_freq/DMA_freq;

SRRC = rcosine(1,M_1,'fir/sqrt',0.3,5);
SRRC = SRRC/max(SRRC);

%% a
x_3;

s3a_demo = x_3.*exp(1i*2*pi*Fc/DAC_freq.*[1:length(x_3)]);

s3a_IIR = filter(p_02_IIR,s3a_demo);

delay = 1;
s3a_down = s3a_IIR([delay:M_2:end]);

s3a_srrc = conv(s3a_down,SRRC,"same");

delay = 5;
a_3ha = s3a_srrc([delay:M_1:end]);

a_3ha = (a_3ha / sqrt(mean(abs(real(a_3ha)).^2)));

rece_angle = angle(a_3ha(6:end-5));
CFO = rece_angle - [0 rece_angle(1:end-1)];

figure;
stem([1:length(a_3ha)],real(a_3ha));

figure;
plot(a_3ha(6:end-5),"O")

%% b
A_3ha = fftshift(fft(real(a_3ha(6:end-5))));
figure;
plot(linspace(-1/2*DAC_freq,1/2*DAC_freq,length(A_3ha)),abs(A_3ha)*DAC_freq);

original = (angle(a_3b));
received = (angle(a_3ha(6:end-5)));
CFO_3b = mean(received-original);

%% c
x_3;

s3a_demo = x_3.*exp(1i*2*pi*(Fc/DAC_freq).*[1:length(x_3)] ).*exp(1i*2*pi*(-CFO_3b));

s3a_IIR = filter(p_02_IIR,s3a_demo);

delay = 1;
s3a_down = s3a_IIR([delay:M_2:end]);

s3a_srrc = conv(s3a_down,SRRC,"same");

delay = 5;
a_3hc = s3a_srrc([delay:M_1:end]);

a_3hc = (a_3hc / sqrt(mean(abs(real(a_3hc)).^2)));

figure;
stem([1:length(a_3hc)],real(a_3hc));

figure;
plot(a_3hc(6:end-5),"O")