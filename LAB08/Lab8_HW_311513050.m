clear all; clc; close all;
%% HW step1
% 1. Let the symbol rate for a system be 1MHz, the sampling
%       rate of the DAC be 4MHz, and sampling rate for DMA filter be 32MHz.

symbol_rate = 1e6;
DAC_rate = 4e6;
DMA_rate = 32e6;
carrier_freq_analog = 8e6;

M_1 = DAC_rate/symbol_rate;
M_2 = DMA_rate/DAC_rate;

n_1 = [-3*M_1:3*M_1] + 1e-6;% Avoid Singularity

%% HW step2
% 2.Let the modulation be BPSK, and the digital pulse shaping is SRRC.

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

%% HW step3
% 3.Design an IIR DMA filter with 5 coefficients that maximizes stopband attenuation.
        
s_03_up_2 = up_sample(M_2,s_02_SRRC_1);       

s_04_IIR = filter(Lab8_hw_IIR,s_03_up_2);

s_04_IIR = filter(Lab8_hw_IIR,s_04_IIR);

S_04_IIR = fftshift(fft(s_04_IIR));



%% HW step4
% 4.Conduct the transmit and receive operation for a sequence (without modulation)
s_05_ch_add_noise = s_04_IIR + randn(1,length(s_04_IIR)).*sqrt(0);

s_06_filter = filter(Lab8_hw_IIR,s_05_ch_add_noise);

s_07_IIR_down_1 = down_sample(M_2,s_06_filter);

s_08_IIR_SRRC1 = conv(s_07_IIR_down_1,SRRC_1_n);
s_08_IIR_SRRC1 = s_08_IIR_SRRC1([floor((length(s_08_IIR_SRRC1)-length(s_07_IIR_down_1))/2)+1 :...
            floor((length(s_08_IIR_SRRC1)-length(s_07_IIR_down_1))/2)+length(s_07_IIR_down_1)]);% aligning

delay = 1;
a_hat_IIR_09 = s_08_IIR_SRRC1([mod(delay,M_1)+1:M_1:length(s_08_IIR_SRRC1)]);

a_hat_IIR_normalized = a_hat_IIR_09/(sqrt(mean(a_hat_IIR_09.^2)));

figure;
stem([1:length(a_hat_IIR_normalized)],real(a_hat_IIR_normalized));
hold on;
stem([1:length(s_BPSK)],real(s_BPSK));
hold off;
legend("a hat IIR, Receive","s BPSK, Transmit");
title_text = "Transmit and Receive operation for a sequence";
title(title_text,"fontsize",12);

%% HW step5
w_c = carrier_freq_analog/symbol_rate/M_1/M_2*2*pi;
carrier = exp(i*w_c.*[1:length(s_04_IIR)]);

s_05_channel = real(s_04_IIR .* carrier);

S_05_channel = fftshift(fft(s_05_channel));

figure;
plot(linspace(-1/2*DMA_rate,1/2*DMA_rate,length(S_04_IIR)),abs(S_04_IIR));
hold on;
plot(linspace(-1/2*DMA_rate,1/2*DMA_rate,length(S_05_channel)),abs(S_05_channel));
hold off;
legend("original","carried");
title_text = "s_c_h Frequency Domain, carrier Freq = "+num2str(carrier_freq_analog);
title(title_text,"fontsize",12);
