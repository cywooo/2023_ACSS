clear all; clc; close all;

symbol_rate = 1e6;
DAC_rate = 16e6;
DMA_rate = 32e6;
carrier_freq_analog = 8e6;
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
s_02_SRRC_1_ = conv(s_01_up_1,SRRC_1_n,"same");
s_02_SRRC_1 = s_02_SRRC_1([floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+1 :...
            floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning
        
s_03_up_2 = up_sample(M_2,s_02_SRRC_1);       

s_04_IIR = filter(Lab9_demo_IIR_2,s_03_up_2);

s_04_IIR = filter(Lab9_demo_IIR_2,s_04_IIR);

S_04_IIR = fftshift(fft(s_04_IIR));

%% practice 01

w_c = carrier_freq_analog/symbol_rate/M_1/M_2*2*pi;

phase_shift = 0.51*pi;
carrier_tran = exp(i*w_c.*[1:length(s_04_IIR)]);
carrier_rece = exp(i*w_c.*[1:length(s_04_IIR)]);
carrier_rece_shift = exp(i*(w_c.*[1:length(s_04_IIR)]+ phase_shift));

s_05_channel = real(s_04_IIR .* carrier_tran);

S_05_channel = fftshift(fft(s_05_channel));

%{
figure;
plot(linspace(-1/2*DMA_rate,1/2*DMA_rate,length(S_04_IIR)),abs(S_04_IIR));
hold on;
plot(linspace(-1/2*DMA_rate,1/2*DMA_rate,length(S_05_channel)),abs(S_05_channel));
hold off;
legend("original","carried");
title_text = "s_c_h Frequency Domain, carrier Freq = "+num2str(carrier_freq_analog);
title(title_text,"fontsize",12);
%}

s_06_rece = s_05_channel ./ carrier_rece;

% for p1
s_07_DMA = filter(Lab9_demo_IIR_2,s_06_rece);

delay = 8;
s_08_reDig = s_07_DMA([mod(delay,M_re_p1)+1:M_re_p1:length(s_07_DMA)]);

s_07_reDig_normalized = s_08_reDig/(sqrt(mean(s_08_reDig.^2)));

figure;
stem([1:length(s_07_reDig_normalized)],real(s_07_reDig_normalized));
hold on;
stem([1:length(s_BPSK)],real(s_BPSK));
hold off;
legend("a hat IIR, Receive","s BPSK, Transmit");
title_text = "Transmit and Receive, practice 01";
title(title_text,"fontsize",12);

%% practice 02
% without shift
s_07_DMA = filter(Lab9_demo_IIR_2,s_06_rece);

delay = 1;
s_08_reDig_p2 = s_07_DMA([mod(delay,M_re_p2)+1:M_re_p2:length(s_07_DMA)]);

s_09_DigFilter = conv(s_08_reDig_p2,SRRC_1_n);
s_09_DigFilter = s_09_DigFilter([floor((length(s_09_DigFilter)-length(s_08_reDig_p2))/2)+1 :...
            floor((length(s_09_DigFilter)-length(s_08_reDig_p2))/2)+length(s_08_reDig_p2)]);% aligning

delay = 2;
s_10_reDig_p2 = s_09_DigFilter([mod(delay,M_re2_p2)+1:M_re2_p2:length(s_09_DigFilter)]);

s_10_reDig_p2_normalized = s_10_reDig_p2/(sqrt(mean(s_10_reDig_p2.^2)));

% with shift
s_06_rece_sh = s_05_channel ./ carrier_rece_shift;

s_07_DMA_sh = filter(Lab9_demo_IIR_2,s_06_rece_sh);

delay = 1;
s_08_reDig_p2_sh = s_07_DMA_sh([mod(delay,M_re_p2)+1:M_re_p2:length(s_07_DMA_sh)]);

s_09_DigFilter_sh = conv(s_08_reDig_p2_sh,SRRC_1_n);
s_09_DigFilter_sh = s_09_DigFilter_sh([floor((length(s_09_DigFilter_sh)-length(s_08_reDig_p2_sh))/2)+1 :...
            floor((length(s_09_DigFilter_sh)-length(s_08_reDig_p2_sh))/2)+length(s_08_reDig_p2_sh)]);% aligning

delay = 2;
s_10_reDig_p2_sh = s_09_DigFilter_sh([mod(delay,M_re2_p2)+1:M_re2_p2:length(s_09_DigFilter_sh)]);

s_10_reDig_p2_sh_normalized = s_10_reDig_p2_sh/(sqrt(mean(s_10_reDig_p2_sh.^2)));

figure;
stem([1:length(s_10_reDig_p2_normalized)],real(s_10_reDig_p2_normalized),"b+");
hold on;
stem([1:length(s_10_reDig_p2_sh_normalized)],real(s_10_reDig_p2_sh_normalized),"kO");
stem([1:length(s_BPSK)],real(s_BPSK),"r*",'LineStyle','none');
hold off;
legend("Receive no shift","Receive with shift","s BPSK, Transmit");
title_text = "Transmit and Receive, practice 02 ";
title(title_text,"fontsize",12);

%=================================
S_06_rece = fftshift(fft(s_06_rece));
S_07_DMA_sh = fftshift(fft(s_07_DMA_sh));
S_08_reDig_p2_sh = fftshift(fft(s_08_reDig_p2_sh));



%=================================
%% practice 03
s_06_IF_rece = s_05_channel .* cos(2*pi*(carrier_freq_analog-IF_freq)/DMA_rate.*[1:length(s_05_channel)]);

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

figure;
stem([1:length(s_11_IF_reDig_normalized)],real(s_11_IF_reDig_normalized),"rO");
hold on;
stem([1:length(s_BPSK)],real(s_BPSK),"b+");
hold off;
legend("IF Receive","s BPSK, Transmit");
title_text = "Transmit and Receive, practice 03 ";
title(title_text,"fontsize",12);
