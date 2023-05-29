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
signal_length = 150;
% signal generate------
s_1 = [sign(rand(1,signal_length)-0.5)];
%sum(1.*(s_1==1))
% up sample by 16-------
s_01_up_1 = up_sample(M_1,s_1);
% Gussian filter ---------
GC_1 = 1;
BT = 0.5;
Gn_1 = [-(M_1-1):(M_1-1)];
%Gussian_filter = GC_1*exp(-1*((2*pi^2)/log10(2))*(BT./ M_1)^2.*Gn_1.^2);
%Gussian_filter = Gussian_filter*(BT)*sqrt(2*pi/log10(2));

Gussian_filter = GC_1*exp(-1*((2*pi^2)/log(2))*(BT./ M_1)^2.*Gn_1.^2);

s_02_Gussianf = conv(s_01_up_1,Gussian_filter);
s_02_Gussianf = s_02_Gussianf([floor((length(s_02_Gussianf)-length(s_01_up_1))/2)+1 :...
            floor((length(s_02_Gussianf)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning

% samation----------------
fd = 150e3;
temp = 0;
for k = [1:length(s_02_Gussianf)]
    temp = s_02_Gussianf(k) + temp;
	sum_s_02_Gussianf(k) = temp;
end

s_03_samation = exp(1i*2*pi*fd*(1/symbol_rate)/M_1*sum_s_02_Gussianf);

% load IF frequency-------
s_04_IF = s_03_samation .* exp(1i*2*pi*IF_freq/DAC_rate.*[1:length(s_03_samation)]);
s_04_IF_real = real(s_04_IF);

% up sample to analog-------
s_05_up_analog = up_sample(M_2,s_04_IF_real);

% quan
% max(s_05_up_analog)=1 min(s_05_up_analog)=-1
nob = 10;
nob = min(nob,12);
DR = 2^0; %+-
code_map = -DR:((DR*2)/(2^nob)):(DR-(DR*2)/(2^nob));
s_05_up_analog = F_point_decision(s_05_up_analog,code_map);
S_05_up_analog = fftshift(fft(s_05_up_analog));

% LPF-----------------------
s_06_toChannel = filter(Lab15_demo_IIR,s_05_up_analog);
S_06_toChannel = fftshift(fft(s_06_toChannel));

[px,f] = pwelch(s_06_toChannel,[],[],[],1);

figure;
plot(f*ANL_rate,10*log10(px));
hold on;
plot(Spectrum_mask_x,Spectrum_mask_y + max(10*log10(px)));
hold off;
title_text = "Transmit Power dB";
title(title_text,"fontsize",12);

figure;
subplot(2,1,1);
plot(linspace(-1/2*ANL_rate,1/2*ANL_rate,length(S_05_up_analog)),abs(S_05_up_analog));
title_text = "up sample";
title(title_text,"fontsize",12);
subplot(2,1,2);
plot(linspace(-1/2*ANL_rate,1/2*ANL_rate,length(S_06_toChannel)),abs(S_06_toChannel));
title_text = "after LPF";
title(title_text,"fontsize",12);

s_07_Channel = s_06_toChannel;
% LPF-----------------------
s_08_receivedLPF = filter(Lab15_demo_IIR,s_07_Channel);
% down sample--------------
s_09_down = s_08_receivedLPF([1:M_2:length(s_08_receivedLPF)]);

% demodulate IF frequency-------
s_10_deIFmodu = s_09_down .* exp(-1i*2*pi*IF_freq/DAC_rate.*[1:length(s_09_down)]);

% SRRC----------------------
%=========
M_srrc = M_1;
n_srrc = [-2*M_srrc:2*M_srrc] + 1e-6;% Avoid Singularity
a = 0.8 ;

A = cos((1+a).*pi.*n_srrc./M_srrc);
B = M_srrc.*sin((1-a).*pi.*n_srrc./M_srrc)./(4.*a.*n_srrc);
C = 1-(4.*a.*n_srrc./M_srrc).^2;
SRRC_n = (4.*a./pi).*(A+B)./C;
%=========
s_11_SRRC = conv(s_10_deIFmodu,SRRC_n);
s_11_SRRC = s_11_SRRC([floor((length(s_11_SRRC)-length(s_10_deIFmodu))/2)+1 :...
            floor((length(s_11_SRRC)-length(s_10_deIFmodu))/2)+length(s_10_deIFmodu)]);% aligning
% take phase------------------
s_12_phase = unwrap(angle(s_11_SRRC))/2/pi/fd*DAC_rate*M_1/2;

figure;
plot([1:length(s_12_phase)],s_12_phase);
title_text = "phase";
title(title_text,"fontsize",12);

% diff-----------------------
s_13_Diff = s_12_phase - [0 s_12_phase(1:end-1)] ;

% Gussian filter ---------
s_14_Gussianf = conv(s_13_Diff,Gussian_filter);
s_14_Gussianf = s_14_Gussianf([floor((length(s_14_Gussianf)-length(s_13_Diff))/2)+1 :...
            floor((length(s_14_Gussianf)-length(s_13_Diff))/2)+length(s_13_Diff)]);% aligning
        
% down sample by 16-------------
delay = 5;
s_15_reDig = s_14_Gussianf([delay:M_1:length(s_14_Gussianf)]);


figure;
stem([1:signal_length-2],(s_15_reDig([2:end-1]))/40,"bo");
hold on;
stem([2:signal_length],s_1([2:end]),"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(real), HW";
title(title_text,"fontsize",12);
