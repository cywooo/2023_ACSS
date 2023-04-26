clear all; clc; close all;

symbol_rate = 1e6;
DAC_rate = 16e6;
DMA_rate = 32e6;
IF_freq = 4e6;

M_1 = DAC_rate/symbol_rate;
M_2 = DMA_rate/DAC_rate;

% n domain SRRC
n_1 = [-5*M_1:5*M_1] + 1e-6;% Avoid Singularity
a = 0.9+1e-6 ;% Avoid Singularity
SRRC_1_n = (4.*a./pi).*(cos((1+a).*pi.*n_1./M_1)+ M_1.*sin((1-a).*pi.*n_1./M_1)./(4.*a.*n_1))...
            ./(1-(4.*a.*n_1./M_1).^2);

imbalance_g = 1;
imbalance_phi = 0;%pi/180*90*-1; %  

%% HW
% transmiter from LAB9
signal_length = 10;

s_real = sign(randn(1,signal_length));
s_imag = sign(randn(1,signal_length));

s_01_up_1_real = up_sample(M_1,s_real); % rate = symbol_rate*M
s_01_up_1_imag = up_sample(M_1,s_imag); % rate = symbol_rate*M

s_02_SRRC_1_real = conv(s_01_up_1_real,SRRC_1_n);
s_02_SRRC_1_real = s_02_SRRC_1_real([floor((length(s_02_SRRC_1_real)-length(s_01_up_1_real))/2)+1 :...
            floor((length(s_02_SRRC_1_real)-length(s_01_up_1_real))/2)+length(s_01_up_1_real)]);% aligning
        
s_02_SRRC_1_imag = conv(s_01_up_1_imag,SRRC_1_n);
s_02_SRRC_1_imag = s_02_SRRC_1_imag([floor((length(s_02_SRRC_1_imag)-length(s_01_up_1_imag))/2)+1 :...
            floor((length(s_02_SRRC_1_imag)-length(s_01_up_1_imag))/2)+length(s_01_up_1_imag)]);% aligning

s_03_up_2_real = up_sample(M_2,s_02_SRRC_1_real);   
s_03_up_2_imag = up_sample(M_2,s_02_SRRC_1_imag);   

s_04_IIR_real = filter(Lab10_demo_IIR,s_03_up_2_real);
%s_04_IIR_real = filter(Lab10_demo_IIR,s_04_IIR_real);
s_04_IIR_imag = filter(Lab10_demo_IIR,s_03_up_2_imag);
%s_04_IIR_imag = filter(Lab10_demo_IIR,s_04_IIR_imag);

xb_real = s_04_IIR_real;
xb_imag = s_04_IIR_imag;

fc = 1/4*DMA_rate; 
x = xb_real.*sqrt(2).*cos(2*pi*fc.*[1:length(xb_real)]/DMA_rate)...
    -xb_imag.*sqrt(2).*sin(2*pi*fc.*[1:length(xb_imag)]/DMA_rate + imbalance_phi).*imbalance_g;

s_06_rece_real = x .* sqrt(2).*cos(2*pi*(fc-IF_freq)/DMA_rate.*[1:length(x)]);
s_07_DMA_real = filter(Lab10_demo_IIR,s_06_rece_real); 
s_07_DMA_real = real(s_07_DMA_real);

s_06_rece_imag = x .* -1.* sqrt(2).*sin(2*pi*(fc-IF_freq)/DMA_rate.*[1:length(x)]);
s_07_DMA_imag = filter(Lab10_demo_IIR,s_06_rece_imag); 
s_07_DMA_imag = real(s_07_DMA_imag);

delay = 1;
s_08_reDig_real = s_07_DMA_real([delay:M_2:length(s_07_DMA_real)]);
s_08_reDig_imag = s_07_DMA_real([delay:M_2:length(s_07_DMA_imag)]);

s_09_IF_toDigFilter_real = s_08_reDig_real .* exp(-1i*2*pi*IF_freq/DAC_rate.*[1:length(s_08_reDig_real)]);
s_09_IF_toDigFilter_imag = s_08_reDig_imag .* exp(-1i*2*pi*IF_freq/DAC_rate.*[1:length(s_08_reDig_imag)]);

s_10_IF_DigFilter_real = conv(s_09_IF_toDigFilter_real,SRRC_1_n);
s_10_IF_DigFilter_real = s_10_IF_DigFilter_real([floor((length(s_10_IF_DigFilter_real)-length(s_09_IF_toDigFilter_real))/2)+1 :...
            floor((length(s_10_IF_DigFilter_real)-length(s_09_IF_toDigFilter_real))/2)+length(s_09_IF_toDigFilter_real)]);% aligning
        
s_10_IF_DigFilter_imag = conv(s_09_IF_toDigFilter_imag,SRRC_1_n);
s_10_IF_DigFilter_imag = s_10_IF_DigFilter_imag([floor((length(s_10_IF_DigFilter_imag)-length(s_09_IF_toDigFilter_imag))/2)+1 :...
            floor((length(s_10_IF_DigFilter_imag)-length(s_09_IF_toDigFilter_imag))/2)+length(s_09_IF_toDigFilter_imag)]);% aligning        
        

delay = 11;
s_11_IF_reDig_real = s_10_IF_DigFilter_real([delay:M_1:length(s_10_IF_DigFilter_real)]);
s_11_IF_reDig_imag = s_10_IF_DigFilter_imag([delay:M_1:length(s_10_IF_DigFilter_imag)]);

s_11_IF_reDig_real_normalized = s_11_IF_reDig_real/(sqrt(mean(abs(s_11_IF_reDig_real).^2)));
s_11_IF_reDig_imag_normalized = s_11_IF_reDig_imag/(sqrt(mean(abs(s_11_IF_reDig_imag).^2)));

figure(1);
subplot(2,1,1);
stem([1:length(s_11_IF_reDig_real_normalized)],real(s_11_IF_reDig_real_normalized),"bo");
hold on;
stem([1:length(s_real)],s_real,"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(real), HW g = "+num2str(imbalance_g)+", phi = "+num2str(imbalance_phi);
title(title_text,"fontsize",12);
subplot(2,1,2);
stem([1:length(s_11_IF_reDig_imag_normalized)],imag(s_11_IF_reDig_imag_normalized),"bo");
hold on;
stem([1:length(s_imag)],s_imag,"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(imag), HW g = "+num2str(imbalance_g)+", phi = "+num2str(imbalance_phi);
title(title_text,"fontsize",12);

