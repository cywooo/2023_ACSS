clear all; clc; close all;

symbol_rate = 1e6;
DAC_rate = 16e6;
DMA_rate = 32e6;

M_1 = DAC_rate/symbol_rate;
M_2 = DMA_rate/DAC_rate;

M_Direct = DMA_rate/symbol_rate;

n_1 = [-2*M_1:2*M_1] + 1e-6;% Avoid Singularity

a = 0.95+1e-6 ;% Avoid Singularity

% n domain SRRC
SRRC_1_n = (4.*a./pi).*(cos((1+a).*pi.*n_1./M_1)+ M_1.*sin((1-a).*pi.*n_1./M_1)./(4.*a.*n_1))...
            ./(1-(4.*a.*n_1./M_1).^2);
        

imbalance_g = 1.5;
imbalance_phi = pi/180*15; %  0;%

%% pratice 01
% transmiter from LAB9
signal_length = 30;
s_1 = sign(randn(1,signal_length))+1i*sign(randn(1,signal_length));

%{
n = [1:201];
x_I = 0 + 1 .*(n >= 50).*(n < 151);
x_Q = 0 + (n-50)/50 .*(n >= 50).*(n < 101)+(151-n)/50.*(n >= 101).*(n < 151);
s_1 = x_I + x_Q*1i;
%}

% upsampling factor is 4
s_01_up_1 = up_sample(M_1,s_1); % rate = symbol_rate*M

% SRRC
s_02_SRRC_1 = conv(s_01_up_1,SRRC_1_n);
s_02_SRRC_1 = s_02_SRRC_1([floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+1 :...
            floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning
        
s_03_up_2 = up_sample(M_2,s_02_SRRC_1);       

% SRRC

s_04_IIR = filter(Lab10_demo_IIR,s_03_up_2);
%s_04_IIR = filter(Lab9_demo_IIR_2,s_04_IIR);
%=============================================
% signal formation IQ imbalance lab10 p.3
xb_real = real(s_04_IIR);
xb_imag = imag(s_04_IIR);

fc = 1/4*DMA_rate; 

x_ideal = xb_real.*sqrt(2).*cos(2*pi*fc.*[1:length(xb_real)]/DMA_rate)-...
    xb_imag.*sqrt(2).*sin(2*pi*fc.*[1:length(xb_imag)]/DMA_rate);

x = xb_real.*sqrt(2).*cos(2*pi*fc.*[1:length(xb_real)]/DMA_rate)-...
    imbalance_g.*xb_imag.*sqrt(2).*sin(2*pi*fc.*[1:length(xb_imag)]/DMA_rate + imbalance_phi);
%=============================================
carrier_rece = exp(1i.*2*pi.*fc.*[1:length(x)]/DMA_rate);

s_06_rece = x ./ carrier_rece;
s_07_DMA = filter(Lab10_demo_IIR,s_06_rece);
            

s_06_rece_ideal = x_ideal ./ carrier_rece;
s_07_DMA_ideal = filter(Lab10_demo_IIR,s_06_rece_ideal);      

theo_ = (1/2*(1+imbalance_g.*exp(1i*imbalance_phi))).*s_04_IIR + (1/2*(1-imbalance_g.*exp(1i*imbalance_phi))).*conj(s_04_IIR);
theo_ = theo_/sum(abs(theo_))*sum(abs(s_07_DMA));

figure;
subplot(2,1,1);
plot([1:length(s_07_DMA)],real(s_07_DMA),"b");
hold on;
plot([1:length(s_07_DMA_ideal)],real(s_07_DMA_ideal),"r");
plot([1:length(theo_)],real(theo_),"k--");
hold off;
legend("imbalance","ideal","theo");
title_text = "channel received compare (real), practice 1 g = "+num2str(imbalance_g)+", phi = "+num2str(imbalance_phi);
title(title_text,"fontsize",12);

subplot(2,1,2);
plot([1:length(s_07_DMA)],imag(s_07_DMA),"b");
hold on;
plot([1:length(s_07_DMA_ideal)],imag(s_07_DMA_ideal),"r");
plot([1:length(theo_)],imag(theo_),"k--");
hold off;
legend("imbalance","ideal","theo");
title_text = "channel received compare (imag), practice 1 g = "+num2str(imbalance_g)+", phi = "+num2str(imbalance_phi);
title(title_text,"fontsize",12);

delay = 11;
s_08_reDig_p1 = s_07_DMA([mod(delay,M_Direct)+1:M_Direct:length(s_07_DMA)]);

s_08_reDig_p1_normalized = real(s_08_reDig_p1)/(sqrt(mean(real(s_08_reDig_p1).^2)))+...
    imag(s_08_reDig_p1)/(sqrt(mean(imag(s_08_reDig_p1).^2))).*1i;

figure;
subplot(2,1,1);
stem([1:length(s_08_reDig_p1_normalized)],real(s_08_reDig_p1_normalized),"bo");
hold on;
stem([1:length(s_1)],real(s_1),"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(real), practice 1 g = "+num2str(imbalance_g)+", phi = "+num2str(imbalance_phi);
title(title_text,"fontsize",12);
subplot(2,1,2);
stem([1:length(s_08_reDig_p1_normalized)],imag(s_08_reDig_p1_normalized),"bo");
hold on;
stem([1:length(s_1)],imag(s_1),"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(imag), practice 1 g = "+num2str(imbalance_g)+", phi = "+num2str(imbalance_phi);
title(title_text,"fontsize",12);

%% pratice 02
% transmiter from LAB9
signal_length = 30;
s_1 = sign(randn(1,signal_length))+1i*sign(randn(1,signal_length));

H = [1 -imbalance_g*sin(imbalance_phi) ; 0  imbalance_g*cos(imbalance_phi)];

s_1_b = inv(H) * [real(s_1);imag(s_1)];
s_1_compensate = s_1_b(1,:)+s_1_b(2,:)*1i;
 
% upsampling factor is 4
s_01_up_1 = up_sample(M_1,s_1_compensate); % rate = symbol_rate*M

% SRRC
s_02_SRRC_1 = conv(s_01_up_1,SRRC_1_n);
s_02_SRRC_1 = s_02_SRRC_1([floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+1 :...
            floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning
        
s_03_up_2 = up_sample(M_2,s_02_SRRC_1);       

s_04_IIR = filter(Lab10_demo_IIR,s_03_up_2);
%s_04_IIR = filter(Lab9_demo_IIR_2,s_04_IIR);
%=============================================
% signal formation IQ imbalance lab10 p.3
xb_real = real(s_04_IIR);
xb_imag = imag(s_04_IIR);

fc = 1/4*DMA_rate; 

x_ideal = xb_real.*sqrt(2).*cos(2*pi*fc.*[1:length(xb_real)]/DMA_rate)-...
    xb_imag.*sqrt(2).*sin(2*pi*fc.*[1:length(xb_imag)]/DMA_rate);

x = xb_real.*sqrt(2).*cos(2*pi*fc.*[1:length(xb_real)]/DMA_rate)-...
    imbalance_g.*xb_imag.*sqrt(2).*sin(2*pi*fc.*[1:length(xb_imag)]/DMA_rate + imbalance_phi);
%=============================================
carrier_rece = exp(1i.*2*pi.*fc.*[1:length(x)]/DMA_rate);

s_06_rece = x ./ carrier_rece;
s_07_DMA = filter(Lab10_demo_IIR,s_06_rece);

s_06_rece_ideal = x_ideal ./ carrier_rece;
s_07_DMA_ideal= filter(Lab10_demo_IIR,s_06_rece_ideal);

s_07_DMA_ideal_ = real(s_07_DMA_ideal)/sum(abs(real(s_07_DMA_ideal)))*sum(abs(real(s_07_DMA)))+...
    imag(s_07_DMA_ideal)/sum(abs(imag(s_07_DMA_ideal)))*sum(abs(imag(s_07_DMA))).*1i;

figure;
subplot(2,1,1);
plot([1:length(s_07_DMA)],real(s_07_DMA),"b");
hold on;
plot([1:length(s_07_DMA_ideal_)],real(s_07_DMA_ideal_),"r--");
hold off;
legend("compensate","ideal");
title_text = "compensate received compare (real), practice 2 g = "+num2str(imbalance_g)+", phi = "+num2str(imbalance_phi);
title(title_text,"fontsize",12);
subplot(2,1,2);
plot([1:length(s_07_DMA)],imag(s_07_DMA),"b");
hold on;
plot([1:length(s_07_DMA_ideal_)],imag(s_07_DMA_ideal_),"r--");
hold off;
legend("compensate","ideal");
title_text = "compensate received compare (imag), practice 2 g = "+num2str(imbalance_g)+", phi = "+num2str(imbalance_phi);
title(title_text,"fontsize",12);

delay = 11;
s_08_reDig_p1 = s_07_DMA([mod(delay,M_Direct)+1:M_Direct:length(s_07_DMA)]);

s_08_reDig_p1_normalized = real(s_08_reDig_p1)/(sqrt(mean(real(s_08_reDig_p1).^2)))+...
    imag(s_08_reDig_p1)/(sqrt(mean(imag(s_08_reDig_p1).^2))).*1i;

figure;
subplot(2,1,1);
stem([1:length(s_08_reDig_p1_normalized)],real(s_08_reDig_p1_normalized),"bo");
hold on;
stem([1:length(s_1)],real(s_1),"r+");
hold off;
legend("re","tr");
title_text = "compensate Transmit & Receive(real), practice 2 g = "+num2str(imbalance_g)+", phi = "+num2str(imbalance_phi);
title(title_text,"fontsize",12);
subplot(2,1,2);
stem([1:length(s_08_reDig_p1_normalized)],imag(s_08_reDig_p1_normalized),"bo");
hold on;
stem([1:length(s_1)],imag(s_1),"r+");
hold off;
legend("re","tr");
title_text = "compensate Transmit & Receive(imag), practice 2 g = "+num2str(imbalance_g)+", phi = "+num2str(imbalance_phi);
title(title_text,"fontsize",12);


