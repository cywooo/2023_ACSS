clear all; clc; close all;

symbol_rate = 1e6;
DAC_rate = 16e6;
ANL_rate = 64e6;
IF_freq = 2e6;
carrier_freq_analog = 8e6;

Spectrum_mask_x = linspace(0,ANL_rate/2,1000);
Spectrum_mask_y = -60.*(Spectrum_mask_x<(IF_freq-2.5e6))+... 
                  -50.*(Spectrum_mask_x>=(IF_freq-2.5e6)).*((Spectrum_mask_x<IF_freq-1.5e6))+...
                  -26.*(Spectrum_mask_x>=(IF_freq-1.5e6)).*((Spectrum_mask_x<IF_freq-1e6))+...
                  0.*(Spectrum_mask_x>=(IF_freq-1e6)).*((Spectrum_mask_x<IF_freq+1e6))+...
                  -26.*(Spectrum_mask_x>=(IF_freq+1e6)).*((Spectrum_mask_x<IF_freq+1.5e6))+...
                  -50.*(Spectrum_mask_x>=(IF_freq+1.5e6)).*((Spectrum_mask_x<IF_freq+2.5e6))+...
                  -60.*(Spectrum_mask_x>=(IF_freq+2.5e6));
              
M_1 = DAC_rate/symbol_rate;
M_2 = ANL_rate/DAC_rate;

%% demo
signal_length = 1e3;
% signal generate------
s_1 = sign(rand(1,signal_length)-0.5)+sign(rand(1,signal_length)-0.5)*1i;

s_phase = angle(s_1);

% samation----------------
temp = 0;
for k = [1:length(s_phase)]
    temp = s_phase(k) + temp;
	sum_s_phase(k) = temp;
end

s_I = sqrt(2)*cos(sum_s_phase);
s_Q = sqrt(2)*sin(sum_s_phase);

s_IQ = s_I + s_Q*1i;

% up sample by 16-------
s_01_up_1 = up_sample(M_1,s_IQ);

% SRRC----------------------
%=========
M_srrc = M_1;
n_srrc = [-2*M_srrc:2*M_srrc] + 1e-6;% Avoid Singularity
a = 0.4 ;

A = cos((1+a).*pi.*n_srrc./M_srrc);
B = M_srrc.*sin((1-a).*pi.*n_srrc./M_srrc)./(4.*a.*n_srrc);
C = 1-(4.*a.*n_srrc./M_srrc).^2;
SRRC_n = (4.*a./pi).*(A+B)./C;
%=========
s_02_SRRC = conv(s_01_up_1,SRRC_n);
s_02_SRRC = s_02_SRRC([floor((length(s_02_SRRC)-length(s_01_up_1))/2)+1 :...
            floor((length(s_02_SRRC)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning

% load IF frequency-------
s_03_IF = s_02_SRRC .* exp(1i*2*pi*IF_freq/DAC_rate.*[1:length(s_02_SRRC)]);
s_04_real = real(s_03_IF);

% up sample to analog-------
s_05_up_analog = up_sample(M_2,s_04_real);

SNR_dB = -2;
SNR = 10^(SNR_dB/10);
signal_power = mean(abs(s_05_up_analog).^2);
noise_power = signal_power/SNR;
s_06_add_noise = s_05_up_analog + randn(1,length(s_05_up_analog))*sqrt(noise_power);

% LPF-----------------------
s_06_toChannel = filter(Lab15_HW_p2_IIR,s_06_add_noise);

s_06_channel = s_06_toChannel;

[px,f] = pwelch(s_06_channel,[],[],[],1);

figure;
plot(f*ANL_rate,10*log10(px));
hold on;
plot(Spectrum_mask_x,Spectrum_mask_y + max(10*log10(px)));
hold off;
title_text = "Transmit Power dB";
title(title_text,"fontsize",12);

% down sample by 4 --------------
s_07_down = s_06_channel([1:M_2:length(s_06_channel)]);

% demodulate IF frequency-------
s_08_deIFmodu = s_07_down .* exp(-1i*2*pi*IF_freq/DAC_rate.*[1:length(s_07_down)]);

s_09_SRRC = conv(s_08_deIFmodu,SRRC_n);
s_09_SRRC = s_09_SRRC([floor((length(s_09_SRRC)-length(s_08_deIFmodu))/2)+1 :...
            floor((length(s_09_SRRC)-length(s_08_deIFmodu))/2)+length(s_08_deIFmodu)]);% aligning
        
% down sample by 16-------------
delay = 5;%5
s_10_down = s_09_SRRC([delay:M_1:length(s_09_SRRC)]);

s_11_phase = angle(s_10_down);

% diff-----------------------
s_12_Diff = s_11_phase - [0 s_11_phase(1:end-1)] ;

s_13_QPSK_hat = sqrt(2)*cos(s_12_Diff) + sqrt(2)*sin(s_12_Diff)*1i;

s_13_hat_I = real(s_13_QPSK_hat);
s_13_hat_Q = imag(s_13_QPSK_hat);
s_1_I = real(s_1);
s_1_Q = imag(s_1);

BER = mean(abs(sign(s_13_hat_I([12:end-9]))-s_1_I([11:end-10]))/2)/2 +...
      mean(abs(sign(s_13_hat_Q([12:end-9]))-s_1_Q([11:end-10]))/2)/2;

EVM_dB = 10*log10( mean(abs(s_13_QPSK_hat([12:end-9])-s_1([11:end-10])).^2) / mean(abs(s_1([11:end-10])).^2));
EVM_percent = sqrt(mean(abs(s_13_QPSK_hat([12:end-9])-s_1([11:end-10])).^2) / mean(abs(s_1([11:end-10])).^2));

figure;
scatter(real(s_13_QPSK_hat([11:end-10])),imag(s_13_QPSK_hat([11:end-10])),"bo");
hold on;
scatter(real(s_1([11:end-10])),imag(s_1([11:end-10])),"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive, EVM_d_B= "+num2str(EVM_dB)+", EVM_%= "+num2str(EVM_percent*100)+"%";
title(title_text,"fontsize",12);

figure;
stem([1:50],s_13_hat_I([12:61]),"bo");
hold on;
stem([1:50],s_1_I([11:60]),"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive, BER = "+num2str(BER);
title(title_text,"fontsize",12);

% figure;
% stem([1:length(s_13_hat_I)],s_13_hat_I,"bo");
% hold on;
% stem([2:length(s_1_I)+1],s_1_I,"r+");
% hold off;
% legend("re","tr");
% title_text = "Transmit and Receive, BER = "+num2str(BER);
% title(title_text,"fontsize",12);
