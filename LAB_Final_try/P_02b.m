clear all; clc; close all;
load('data_given');

signal_freq = 1e6;
DMA_freq = 16e6;
DAC_freq = 128e6;
Fc = 32e6;
Fif = 2e6;

SRRC = rcosine(1,16,'fir/sqrt',0.3,5);
SRRC = SRRC/max(SRRC);
M_1 = DMA_freq/signal_freq;
M_2 = DAC_freq/DMA_freq;

guard = 10;
s_01 = [zeros(1,guard) a_2a zeros(1,guard)];
% up 16
s_02_up = zeros(M_1,length(s_01));
s_02_up(1,:) = s_01;
s_02_up = reshape(s_02_up,1,length(s_01)*M_1);

s_03_srrc = conv(s_02_up,SRRC,"same");

% up 8
s_04_up = zeros(M_2,length(s_03_srrc));
s_04_up(1,:) = s_03_srrc;
s_04_up = reshape(s_04_up,1,length(s_03_srrc)*M_2);

s_05_IIR = filter(p_02_IIR,s_04_up);

s_06_carrier = s_05_IIR.*exp(1i*2*pi*(Fc/DAC_freq).*[1:length(s_05_IIR)]);

s_07_real = real(s_06_carrier);

s_08_channel = s_07_real + [zeros(1,guard*M_1*M_2) intf_2b(1:length(a_2a)*M_1*M_2) zeros(1,guard*M_1*M_2)];

s_08_Bandpass = filter(p_02b_IIR_bandpass,s_08_channel);

S_07_real = fftshift(fft(s_07_real));
S_08_channel = fftshift(fft(s_08_channel));
S_08_Bandpass = fftshift(fft(s_08_Bandpass));
figure;
plot(linspace(-1/2*DAC_freq,1/2*DAC_freq,length(S_07_real)),abs(S_07_real)*DAC_freq);
figure;
plot(linspace(-1/2,1/2,length(S_08_channel)),abs(S_08_channel));
figure;
plot(linspace(-1/2,1/2,length(S_08_Bandpass)),abs(S_08_Bandpass));

s_09_cos = s_08_Bandpass.*cos(2*pi*(Fc-Fif)/DAC_freq.*[1:length(s_08_Bandpass)]);

s_10_IIR = filter(p_02_IIR,s_09_cos);

delay = 3;
s_11_down = s_10_IIR([delay:M_2:end]);

s_12_deif = s_11_down.*exp(-1i*2*pi*Fif/DMA_freq.*[1:length(s_11_down)]);

s_13_srrc = conv(s_12_deif,SRRC,"same");

delay = 14;
a_2hb = s_13_srrc([delay:M_1:end]);
a_2hb = a_2hb / sqrt(mean(abs(real(a_2hb([guard+1:end-guard]))).^2));


figure;
stem([1:length(a_2a)],a_2a);
hold on;
stem([1:length(a_2a)],real(a_2hb([guard+1:end-guard])));
hold off;

noise = real(a_2hb([guard+1:end-guard])) - a_2a;
noise_power = mean(abs(noise).^2);
SNR_2b = 10*log10(1/noise_power);
