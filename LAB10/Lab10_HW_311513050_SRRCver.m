clear all; clc; close all;

symbol_rate = 1e6;
DAC_rate = 16e6;
DMA_rate = 32e6;
IF_freq = 4e6;

M_1 = DAC_rate/symbol_rate;
M_2 = DMA_rate/DAC_rate;

% n domain SRRC
n_1 = [-2*M_1:2*M_1] + 1e-6;% Avoid Singularity
n_2 = [-2*M_2:2*M_2] + 1e-6;% Avoid Singularity
a = 0.5 ;% Avoid Singularity

A = cos((1+a).*pi.*n_1./M_1);
B = M_1.*sin((1-a).*pi.*n_1./M_1)./(4.*a.*n_1);
C = 1-(4.*a.*n_1./M_1).^2;
SRRC_1_n = (4.*a./pi).*(A+B)./C;

A = cos((1+a).*pi.*n_2./M_2);
B = M_2.*sin((1-a).*pi.*n_2./M_2)./(4.*a.*n_2);
C = 1-(4.*a.*n_2./M_2).^2;
SRRC_2_n = (4.*a./pi).*(A+B)./C;

imbalance_g = 1;
imbalance_phi = 0;%  2*pi/360*30; %

%% HW
signal_length = 100;
% signal generate------
s_real = sign(randn(1,signal_length)); %[-1 ones(1,signal_length/2-1) -1*ones(1,signal_length/2-1) 1]; %
s_imag = sign(randn(1,signal_length)); %[-1*ones(1,signal_length/2) ones(1,signal_length/2)]; %
s_1 = s_real + s_imag.*1i;
% up sample by 16-------
s_01_up_1_real = up_sample(M_1,s_real); 
s_01_up_1_imag = up_sample(M_1,s_imag);
s_01_up_1 = up_sample(M_1,s_1);
% SRRC filter ---------
s_02_SRRC_1_real = conv(s_01_up_1_real,SRRC_1_n);
s_02_SRRC_1_real = s_02_SRRC_1_real([floor((length(s_02_SRRC_1_real)-length(s_01_up_1_real))/2)+1 :...
            floor((length(s_02_SRRC_1_real)-length(s_01_up_1_real))/2)+length(s_01_up_1_real)]);% aligning
        
s_02_SRRC_1_imag = conv(s_01_up_1_imag,SRRC_1_n);
s_02_SRRC_1_imag = s_02_SRRC_1_imag([floor((length(s_02_SRRC_1_imag)-length(s_01_up_1_imag))/2)+1 :...
            floor((length(s_02_SRRC_1_imag)-length(s_01_up_1_imag))/2)+length(s_01_up_1_imag)]);% aligning
        
s_02_SRRC_1 = conv(s_01_up_1,SRRC_1_n);
s_02_SRRC_1 = s_02_SRRC_1([floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+1 :...
            floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning
        
%s_02_SRRC_1_real = filter(Lab10_HW_IIR_16,s_01_up_1_real);         
%s_02_SRRC_1_imag = filter(Lab10_HW_IIR_16,s_01_up_1_imag);                 
%-----------------------
% up sample by 2-------
s_03_up_2_real = up_sample(M_2,s_02_SRRC_1_real);   
s_03_up_2_imag = up_sample(M_2,s_02_SRRC_1_imag);   
s_03_up_2 = up_sample(M_2,s_02_SRRC_1);   
% DMA filter------------
%s_04_IIR_real = filter(Lab10_demo_IIR,s_03_up_2_real);
%s_04_IIR_imag = filter(Lab10_demo_IIR,s_03_up_2_imag);

s_04_IIR_real = conv(s_03_up_2_real,SRRC_2_n);
s_04_IIR_real = s_04_IIR_real([floor((length(s_04_IIR_real)-length(s_03_up_2_real))/2)+1 :...
            floor((length(s_04_IIR_real)-length(s_03_up_2_real))/2)+length(s_03_up_2_real)]);% aligning
        
s_04_IIR_imag = conv(s_03_up_2_imag,SRRC_2_n);
s_04_IIR_imag = s_04_IIR_imag([floor((length(s_04_IIR_imag)-length(s_03_up_2_imag))/2)+1 :...
            floor((length(s_04_IIR_imag)-length(s_03_up_2_imag))/2)+length(s_03_up_2_imag)]);% aligning
        
s_04_IIR = conv(s_03_up_2,SRRC_2_n);
s_04_IIR = s_04_IIR([floor((length(s_04_IIR)-length(s_03_up_2))/2)+1 :...
            floor((length(s_04_IIR)-length(s_03_up_2))/2)+length(s_03_up_2)]);% aligning

%-----------------------
% IQ imbalance----------
xb_real = s_04_IIR_real;
xb_imag = s_04_IIR_imag;

xb_real = real(s_04_IIR);
xb_imag = imag(s_04_IIR);

fc = 1/4*DMA_rate;
x =  xb_real.*sqrt(2).*cos(2*pi*fc/DMA_rate.*[1:length(xb_real)])...
    -xb_imag.*sqrt(2).*sin(2*pi*fc/DMA_rate.*[1:length(xb_imag)] + imbalance_phi).*imbalance_g;
% IF downconvert to IF freq.----
wc = 2*pi*fc/DMA_rate;
wif = 2*pi*IF_freq/DMA_rate;
s_06_rece = x .* cos( (wc - wif).*[1:length(x)]);
% received DMA------------------
%s_07_DMA = filter(Lab10_demo_IIR,s_06_rece); 

s_07_DMA = conv(s_06_rece,SRRC_2_n);
s_07_DMA = s_07_DMA([floor((length(s_07_DMA)-length(s_06_rece))/2)+1 :...
            floor((length(s_07_DMA)-length(s_06_rece))/2)+length(s_06_rece)]);% aligning

% down sample by 2--------------
s_08_reDig = s_07_DMA([1:M_2:length(s_07_DMA)]);
% downconvert to base freq.-----
wif = 2*pi*IF_freq/DAC_rate;
s_09_IF_toDigFilter = s_08_reDig .* exp(-1i * wif .* [1:length(s_08_reDig)]);
% SRRC filter ------------------
s_10_IF_DigFilter = conv(s_09_IF_toDigFilter,SRRC_1_n);
s_10_IF_DigFilter = s_10_IF_DigFilter([floor((length(s_10_IF_DigFilter)-length(s_09_IF_toDigFilter))/2)+1 :...
            floor((length(s_10_IF_DigFilter)-length(s_09_IF_toDigFilter))/2)+length(s_09_IF_toDigFilter)]);% aligning
%s_10_IF_DigFilter = filter(Lab10_HW_IIR_16,s_09_IF_toDigFilter);         
%-------------------------------
% down sample by 16-------------
delay = 1;
s_11_IF_reDig = s_10_IF_DigFilter([delay:M_1:length(s_10_IF_DigFilter)]);
%-------------------------------
% normalized-------------
s_11_IF_reDig = s_11_IF_reDig* exp(-1i*2*pi*-45/360); %-45
s_11_IF_reDig_normalized = real(s_11_IF_reDig)/(sqrt(mean(real(s_11_IF_reDig).^2)))+...
    imag(s_11_IF_reDig)/(sqrt(mean(imag(s_11_IF_reDig).^2)))*1i;


figure(1);
subplot(2,1,1);
stem([1:signal_length],real(s_11_IF_reDig_normalized),"bo");
hold on;
stem([1:signal_length],s_real,"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(real), HW g = "+num2str(imbalance_g)+", phi = "+num2str(imbalance_phi);
title(title_text,"fontsize",12);

subplot(2,1,2);
stem([1:signal_length],imag(s_11_IF_reDig_normalized),"bo");
hold on;
stem([1:signal_length],s_imag,"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(imag), HW g = "+num2str(imbalance_g)+", phi = "+num2str(imbalance_phi);
title(title_text,"fontsize",12);

figure;
scatter(real(s_11_IF_reDig_normalized),imag(s_11_IF_reDig_normalized),"b");
hold on;
scatter(s_real,s_imag,"r");
%scatter(s_real(10),s_imag(10),"+R");
%scatter(real(s_11_IF_reDig_normalized(10)),imag(s_11_IF_reDig_normalized(10)),"XK");
%scatter(s_real(20),s_imag(20),"+K");
%scatter(real(s_11_IF_reDig_normalized(20)),imag(s_11_IF_reDig_normalized(20)),"XR");
hold off;

figure;
scatter(real(s_02_SRRC_1),imag(s_02_SRRC_1),"b");
title_text = "1";
title(title_text,"fontsize",12);

figure;
scatter(real(s_04_IIR),imag(s_04_IIR),"b");
title_text = "2";
title(title_text,"fontsize",12);

figure;
scatter(real(s_07_DMA),imag(s_07_DMA),"b");
title_text = "3";
title(title_text,"fontsize",12);

figure;
scatter(real(s_10_IF_DigFilter),imag(s_10_IF_DigFilter),"b");
title_text = "4";
title(title_text,"fontsize",12);

figure;
scatter(real(s_09_IF_toDigFilter),imag(s_09_IF_toDigFilter),"b");
title_text = "pass exp";
title(title_text,"fontsize",12);