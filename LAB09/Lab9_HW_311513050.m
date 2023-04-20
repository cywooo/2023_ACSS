clear all; clc; close all;

symbol_rate = 1e6;
DAC_rate = 16e6;
DMA_rate = 64e6;
carrier_freq_analog = 16e6;
IF_freq = 4e6;

ADC_rate_p1 = 1e6;
ADC_rate_p2 = 16e6;
ADC_rate_p3 = 16e6;


M_1 = DAC_rate/symbol_rate;
M_2 = DMA_rate/DAC_rate;

M_re_p1 = DMA_rate/ADC_rate_p1;

M_re_p2 = DMA_rate/ADC_rate_p2;
M_re2_p2 = ADC_rate_p2/symbol_rate;

n_1 = [-3*M_1:3*M_1] + 1e-6;% Avoid Singularity

a = 0.5+1e-6 ;% Avoid Singularity
% n domain SRRC
SRRC_1_n = (4.*a./pi).*(cos((1+a).*pi.*n_1./M_1)+ M_1.*sin((1-a).*pi.*n_1./M_1)./(4.*a.*n_1))...
            ./(1-(4.*a.*n_1./M_1).^2);

bpsk_length = 1e2;
s_BPSK = sign(randn(1,bpsk_length));

% upsampling factor is 4
s_01_up_1 = up_sample(M_1,s_BPSK); % rate = symbol_rate*M

% SRRC
s_02_SRRC_1 = conv(s_01_up_1,SRRC_1_n);
s_02_SRRC_1 = s_02_SRRC_1([floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+1 :...
            floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning
        
s_03_up_2 = up_sample(M_2,s_02_SRRC_1);       

s_04_IIR = filter(Lab9_demo_IIR_2,s_03_up_2);

s_04_IIR = filter(Lab9_demo_IIR_2,s_04_IIR);

S_04_IIR = fftshift(fft(s_04_IIR));

%% practice 01

w_c = carrier_freq_analog/symbol_rate/M_1/M_2*2*pi;

phase_shift = 0*pi;
carrier_tran = exp(i*w_c.*[1:length(s_04_IIR)]);
carrier_rece = exp(i*w_c.*[1:length(s_04_IIR)]);
carrier_rece_shift = exp(i*(w_c.*[1:length(s_04_IIR)]+ phase_shift));

s_05_channel = real(s_04_IIR .* carrier_tran);

s_05_channel_noise = s_05_channel + randn(1,length(s_05_channel))*0.5;


s_06_IF_rece = s_05_channel_noise .* cos(2*pi*(carrier_freq_analog-IF_freq)/DMA_rate.*[1:length(s_05_channel_noise)]);

s_07_IF_DMA = filter(Lab9_demo_IIR_2,s_06_IF_rece);

delay = 1;
s_08_IF_reDig = s_07_IF_DMA([mod(delay,M_re_p2)+1:M_re_p2:length(s_07_IF_DMA)]);

w_if = IF_freq/DAC_rate*2*pi;
s_09_IF_toDigFilter = s_08_IF_reDig .* exp(-i*w_if.*[1:length(s_08_IF_reDig)]);

s_10_IF_DigFilter = conv(s_09_IF_toDigFilter,SRRC_1_n);
s_10_IF_DigFilter = s_10_IF_DigFilter([floor((length(s_10_IF_DigFilter)-length(s_09_IF_toDigFilter))/2)+1 :...
            floor((length(s_10_IF_DigFilter)-length(s_09_IF_toDigFilter))/2)+length(s_09_IF_toDigFilter)]);% aligning

delay = 2;
s_11_IF_reDig = s_10_IF_DigFilter([mod(delay,M_re2_p2)+1:M_re2_p2:length(s_10_IF_DigFilter)]);

s_11_IF_reDig_normalized = s_11_IF_reDig/(sqrt(mean(s_11_IF_reDig.^2)));
s_11_IF_reDig_detect = sign(real(s_11_IF_reDig))+sign(imag(s_11_IF_reDig))*1i;
figure;
stem([1:length(s_11_IF_reDig_detect)],real(s_11_IF_reDig_detect),"rO");
hold on;
stem([1:length(s_BPSK)],real(s_BPSK),"b+");
hold off;
legend("IF Receive","s BPSK, Transmit");
title_text = "Transmit and Receive";
title(title_text,"fontsize",12);