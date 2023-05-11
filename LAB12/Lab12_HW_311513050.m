clear all; clc; close all;

symbol_rate = 1e6;
DAC_rate = 16e6;
ANL_rate = 64e6;
IF_freq = 2e6;


M_1 = DAC_rate/symbol_rate;
M_2 = ANL_rate/DAC_rate;

%% HW
signal_length = 50 +20;
% signal generate------
s_1 = [zeros(1,10) sign(randn(1,signal_length-20)) zeros(1,10)];
% up sample by 16-------
s_01_up_1 = up_sample(M_1,s_1);
% Gussian filter ---------
fd = 150e3; %150e3;

%B = 75e3; %defult
%BT = B/fd; %BT = 0.5;

BT = 0.48;
B = BT*fd;

GC_1 = 1;
Gn_1 = [-floor(M_1*2-1):floor(M_1*2-1)];
Gussian_filter = GC_1*exp(-1*((2*pi^2)/log(2))*(BT./ M_1)^2.*Gn_1.^2);

figure;
stem(Gn_1,Gussian_filter);
title_text = "Gussian filter, n domain, M = "+num2str(M_1)+" ,BT = "+num2str(BT);
title(title_text,"fontsize",12);

s_02_Gussianf = conv(s_01_up_1,Gussian_filter);
s_02_Gussianf = s_02_Gussianf([floor((length(s_02_Gussianf)-length(s_01_up_1))/2)+1 :...
            floor((length(s_02_Gussianf)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning
% samation----------------
%sum_s_02_Gussianf = s_02_Gussianf * triu(ones(length(s_02_Gussianf)));
for k = [1:length(s_02_Gussianf)]
    sum_s_02_Gussianf(k) = sum(s_02_Gussianf([1:k]));
end

s_03_samation = exp(1i*2*pi*fd*(1/symbol_rate)*sum_s_02_Gussianf);

% load IF frequency-------
s_04_IF = s_03_samation .* exp(1i*2*pi*IF_freq/DAC_rate.*[1:length(s_03_samation)]);
s_04_IF_real = real(s_04_IF);

% up sample to analog-------
s_05_up_analog = up_sample(M_2,s_04_IF_real);  

s_05_up_analog_FeqDomain = fftshift(fft(s_05_up_analog));

% LPF-----------------------
s_06_toChannel = filter(Lab12_demo_IIR,s_05_up_analog);

s_06_toChannel_FeqDomain = fftshift(fft(s_06_toChannel));

figure;
subplot(2,1,1);
plot(linspace(-1/2,1/2,length(s_05_up_analog_FeqDomain)),abs(s_05_up_analog_FeqDomain));
title_text = "up sample";
title(title_text,"fontsize",12);
subplot(2,1,2);
plot(linspace(-1/2,1/2,length(s_06_toChannel_FeqDomain)),abs(s_06_toChannel_FeqDomain));
title_text = "after LPF";
title(title_text,"fontsize",12);

% Channel-------------------
s_07_Channel = s_06_toChannel;
% LPF-----------------------
s_08_receivedLPF = filter(Lab12_demo_IIR,s_07_Channel);
% down sample--------------
delay = 2;
s_09_down = s_08_receivedLPF([delay:M_2:length(s_08_receivedLPF)]);

% demodulate IF frequency-------
s_10_deIFmodu = s_09_down .* exp(-1i*2*pi*IF_freq/DAC_rate.*[1:length(s_09_down)]);

% SRRC----------------------
%=========
M_srrc = M_1;
n_srrc = [-3*M_srrc:3*M_srrc] + 1e-6;% Avoid Singularity
a = 0.5 ;% Avoid Singularity

A_ = cos((1+a).*pi.*n_srrc./M_srrc);
B_ = M_srrc.*sin((1-a).*pi.*n_srrc./M_srrc)./(4.*a.*n_srrc);
C_ = 1-(4.*a.*n_srrc./M_srrc).^2;
SRRC_n = (4.*a./pi).*(A_+B_)./C_;
%=========
s_11_SRRC = conv(s_10_deIFmodu,SRRC_n);
s_11_SRRC = s_11_SRRC([floor((length(s_11_SRRC)-length(s_10_deIFmodu))/2)+1 :...
            floor((length(s_11_SRRC)-length(s_10_deIFmodu))/2)+length(s_10_deIFmodu)]);% aligning
% take phase------------------
s_12_phase = unwrap(angle(s_11_SRRC))/2/pi/fd*DAC_rate;

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
delay = 13;
s_15_reDig = s_14_Gussianf([delay:M_1:length(s_14_Gussianf)]);

s_16_compare = sign(s_15_reDig([11:11+signal_length-21]));
s_1_compare = s_1([11:11+signal_length-21]);

for j=[1:16]
    if j==1; figure; end
    delay_ = j;%j*2-1;
    subplot(4,4,j);
    s_15_reDig_ = s_14_Gussianf([delay_:M_1:length(s_14_Gussianf)]);
    stem([1:signal_length],(s_15_reDig_)/25,"bo");
    hold on;
    stem([1:signal_length],s_1,"r+");
    hold off;
    legend("re","tr");
    title_text = "Transmit and Receive(real), delay = "+num2str(delay_);
    title(title_text,"fontsize",12);
    s_16_compare = sign(s_15_reDig_([11:11+signal_length-21]));
    BER(j) = sum(abs(s_16_compare-s_1_compare))/(signal_length-20);
end
BER_min = min(BER);
optimal_delay = find(BER==BER_min)*2-1;

%{
figure;
stem([1:signal_length],(s_15_reDig)/25,"bo");
hold on;
stem([1:signal_length],s_1,"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(real), HW";
title(title_text,"fontsize",12);
%}
