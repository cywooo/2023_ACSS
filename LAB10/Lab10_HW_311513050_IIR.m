clear all; clc; close all;

symbol_rate = 1e6;
DAC_rate = 16e6;
DMA_rate = 32e6;
IF_freq = 4e6;

M_1 = DAC_rate/symbol_rate;
M_2 = DMA_rate/DAC_rate;

% n domain SRRC
n_1 = [-2*M_1:2*M_1] + 1e-6;% Avoid Singularity
a = 0.9+1e-6 ;% Avoid Singularity
SRRC_1_n = (4.*a./pi).*(cos((1+a).*pi.*n_1./M_1)+ M_1.*sin((1-a).*pi.*n_1./M_1)./(4.*a.*n_1))...
            ./(1-(4.*a.*n_1./M_1).^2);

imbalance_g = 1;
imbalance_phi = 0;%pi/180*15; %  

%% HW
% transmiter from LAB9
signal_length = 50;
%s_1 = sign(randn(1,signal_length))+1i*sign(randn(1,signal_length));
%------------------------------
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
s_04_IIR_real = filter(Lab10_demo_IIR,s_04_IIR_real);
s_04_IIR_imag = filter(Lab10_demo_IIR,s_03_up_2_imag);
s_04_IIR_imag = filter(Lab10_demo_IIR,s_04_IIR_imag);
%------------------------------
%{
% upsampling factor is 4
s_01_up_1 = up_sample(M_1,s_1); % rate = symbol_rate*M

% SRRC
s_02_SRRC_1 = conv(s_01_up_1,SRRC_1_n);
s_02_SRRC_1 = s_02_SRRC_1([floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+1 :...
            floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning
        
s_03_up_2 = up_sample(M_2,s_02_SRRC_1);       

s_04_IIR = filter(Lab10_demo_IIR,s_03_up_2);
s_04_IIR = filter(Lab10_demo_IIR,s_04_IIR);
%}
%=============================================
% signal formation IQ imbalance lab10 p.3
xb_real = s_04_IIR_real; %real(s_04_IIR); %
xb_imag = s_04_IIR_imag; %imag(s_04_IIR); %

XX = fftshift(fft(xb_real+1i*xb_imag));
figure;
plot(linspace(-1/2.*DMA_rate,1/2.*DMA_rate,length(XX)),abs(XX));
%hold on;

fc = 1/4*DMA_rate; 
x = xb_real.*sqrt(2).*cos(2*pi*fc.*[1:length(xb_real)]/DMA_rate)-...
    imbalance_g.*xb_imag.*sqrt(2).*sin(2*pi*fc.*[1:length(xb_imag)]/DMA_rate + imbalance_phi);

%=============================================
s_06_rece = x .* cos(2*pi*(fc-IF_freq)/DMA_rate.*[1:length(x)]);

s_07_DMA = filter(Lab10_demo_IIR,s_06_rece);

%dd
X = fftshift(fft(x));
figure;
plot(linspace(-1/2.*DMA_rate,1/2.*DMA_rate,length(X)),abs(X));
S_06_rece = fftshift(fft(s_06_rece));
figure;
plot(linspace(-1/2.*DMA_rate,1/2.*DMA_rate,length(S_06_rece)),abs(S_06_rece));
S_07_DMA = fftshift(fft(s_07_DMA));
figure;
plot(linspace(-1/2.*DMA_rate,1/2.*DMA_rate,length(S_07_DMA)),abs(S_07_DMA));
%dd   

delay = 1;
s_08_reDig = s_07_DMA([delay:M_2:length(s_07_DMA)]);

S_08_reDig = fftshift(fft(s_08_reDig));
figure;
plot(linspace(-1/2.*DAC_rate,1/2.*DAC_rate,length(S_08_reDig)),abs(S_08_reDig));

w_if = 2*pi*IF_freq/DAC_rate;
s_09_IF_toDigFilter = s_08_reDig .* exp(-1i*w_if.*[1:length(s_08_reDig)]);

S_09_IF_toDigFilter = fftshift(fft(s_09_IF_toDigFilter));
figure;
plot(linspace(-1/2.*DAC_rate,1/2.*DAC_rate,length(S_09_IF_toDigFilter)),abs(S_09_IF_toDigFilter));



s_10_IF_DigFilter = conv(s_09_IF_toDigFilter,SRRC_1_n);
s_10_IF_DigFilter = s_10_IF_DigFilter([floor((length(s_10_IF_DigFilter)-length(s_09_IF_toDigFilter))/2)+1 :...
            floor((length(s_10_IF_DigFilter)-length(s_09_IF_toDigFilter))/2)+length(s_09_IF_toDigFilter)]);% aligning
        
s_10_IF_DigFilter_ = conv(s_10_IF_DigFilter,SRRC_1_n);
s_10_IF_DigFilter_ = s_10_IF_DigFilter_([floor((length(s_10_IF_DigFilter_)-length(s_10_IF_DigFilter))/2)+1 :...
            floor((length(s_10_IF_DigFilter_)-length(s_10_IF_DigFilter))/2)+length(s_10_IF_DigFilter)]);% aligning
%s_10_IF_DigFilter = s_10_IF_DigFilter_;
%s_10_IF_DigFilter = filter(Lab10_demo_IIR,s_09_IF_toDigFilter);

s_123=upsample(s_real,16);
figure;
stem([1:length(s_10_IF_DigFilter)],real(s_10_IF_DigFilter.*abs(real(s_123))));
hold on ;
plot([1:length(s_10_IF_DigFilter)],real(s_10_IF_DigFilter));
stem([1:length(s_123)],real(s_123));

S_10_IF_DigFilter = fftshift(fft(s_10_IF_DigFilter));
figure;
plot(linspace(-1/2.*DAC_rate,1/2.*DAC_rate,length(S_10_IF_DigFilter)),abs(S_10_IF_DigFilter)*max(XX)/max(S_10_IF_DigFilter));

delay = 4;
s_11_IF_reDig = s_10_IF_DigFilter([delay:M_1:length(s_10_IF_DigFilter)]);

s_11_IF_reDig_normalized = real(s_11_IF_reDig)/(sqrt(mean(real(s_11_IF_reDig).^2)))+...
    imag(s_11_IF_reDig)/(sqrt(mean(imag(s_11_IF_reDig).^2))).*1i;

figure(1);
subplot(2,1,1);
stem([1:length(s_11_IF_reDig_normalized)],real(s_11_IF_reDig_normalized),"bo");
hold on;
stem([1:length(s_real)],s_real,"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(real), HW g = "+num2str(imbalance_g)+", phi = "+num2str(imbalance_phi);
title(title_text,"fontsize",12);
subplot(2,1,2);
stem([1:length(s_11_IF_reDig_normalized)],imag(s_11_IF_reDig_normalized),"bo");
hold on;
stem([1:length(s_imag)],s_imag,"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(imag), HW g = "+num2str(imbalance_g)+", phi = "+num2str(imbalance_phi);
title(title_text,"fontsize",12);

