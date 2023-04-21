clear all; clc; close all;

symbol_rate = 1e6;
DAC_rate = 16e6;
DMA_rate = 32e6;

M_1 = DAC_rate/symbol_rate;
M_2 = DMA_rate/DAC_rate;

M_re_p1 = DMA_rate/symbol_rate;

n_1 = [-3*M_1:3*M_1] + 1e-6;% Avoid Singularity
n_2 = [-3*M_2:3*M_2] + 1e-6;% Avoid Singularity

a = 0.5+1e-6 ;% Avoid Singularity

% n domain SRRC
SRRC_1_n = (4.*a./pi).*(cos((1+a).*pi.*n_1./M_1)+ M_1.*sin((1-a).*pi.*n_1./M_1)./(4.*a.*n_1))...
            ./(1-(4.*a.*n_1./M_1).^2);
        
SRRC_2_n = (4.*a./pi).*(cos((1+a).*pi.*n_2./M_2)+ M_2.*sin((1-a).*pi.*n_2./M_2)./(4.*a.*n_2))...
            ./(1-(4.*a.*n_2./M_2).^2);

imbalance_g = 1;
imbalance_f = 0;%pi/180*15;

%% pratice 01
% transmiter from LAB9
signal_length = 1e2;
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
s_04_SRRC_2 = conv(s_03_up_2,SRRC_2_n);
s_04_SRRC_2 = s_04_SRRC_2([floor((length(s_04_SRRC_2)-length(s_03_up_2))/2)+1 :...
                floor((length(s_04_SRRC_2)-length(s_03_up_2))/2)+length(s_03_up_2)]);% aligning

%s_04_IIR = filter(Lab9_demo_IIR_2,s_04_IIR);
%=============================================
% signal formation IQ imbalance lab10 p.3
xb_real = real(s_04_SRRC_2);
xb_imag = imag(s_04_SRRC_2);

fc = 1/4*DMA_rate; 

x_ideal = xb_real.*sqrt(2).*cos(2*pi*fc.*[1:length(xb_real)]/DMA_rate)-...
    xb_imag.*sqrt(2).*sin(2*pi*fc.*[1:length(xb_imag)]/DMA_rate);

x = xb_real.*sqrt(2).*cos(2*pi*fc.*[1:length(xb_real)]/DMA_rate)-...
    imbalance_g.*xb_imag.*sqrt(2).*sin(2*pi*fc.*[1:length(xb_imag)]/DMA_rate + imbalance_f);

%=============================================
carrier_rece = exp(1i.*2*pi.*fc.*[1:length(x)]/DMA_rate);

s_06_rece = x ./ carrier_rece;

s_07 = conv(s_06_rece,SRRC_2_n);
s_07 = s_07([floor((length(s_07)-length(s_06_rece))/2)+1 :...
                floor((length(s_07)-length(s_06_rece))/2)+length(s_06_rece)]);% aligning
delay = 1;
s_08 = s_07([mod(delay,M_2)+1:M_2:length(s_07)]);            

s_09 = conv(s_08,SRRC_1_n);
s_09 = s_09([floor((length(s_09)-length(s_08))/2)+1 :...
                floor((length(s_09)-length(s_08))/2)+length(s_08)]);% aligning

delay = 1;
s_10 = s_09([mod(delay,M_1)+1:M_1:length(s_09)]);     

%%%%
s_06_rece_ideal = x_ideal ./ carrier_rece;

s_07_ideal = conv(s_06_rece_ideal,SRRC_2_n);
s_07_ideal = s_07_ideal([floor((length(s_07_ideal)-length(s_06_rece_ideal))/2)+1 :...
                floor((length(s_07_ideal)-length(s_06_rece_ideal))/2)+length(s_06_rece_ideal)]);% aligning
delay = 1;
s_08_ideal = s_07_ideal([mod(delay,M_2)+1:M_2:length(s_07_ideal)]);            

s_09_ideal = conv(s_08_ideal,SRRC_1_n);
s_09_ideal = s_09_ideal([floor((length(s_09_ideal)-length(s_08_ideal))/2)+1 :...
                floor((length(s_09_ideal)-length(s_08_ideal))/2)+length(s_08_ideal)]);% aligning

delay = 1;
s_10_ideal = s_09_ideal([mod(delay,M_2)+1:M_1:length(s_09_ideal)]);                 
            
%{
theo_ = (1/2*(1+imbalance_g.*exp(1i*imbalance_f))).*s_04_IIR + (1/2*(1-imbalance_g.*exp(1i*imbalance_f))).*conj(s_04_IIR);
theo_ = theo_/sum(abs(theo_))*sum(abs(s_07_DMA));

figure;
plot([1:length(s_07)],real(s_07),"b");
hold on;
plot([1:length(s_07_ideal)],real(s_07_ideal),"r");
plot([1:length(theo_)],real(theo_),"k--");
hold off;
legend("imbalance","ideal","theo");
title_text = "channel received compare (real), practice 01  ";
title(title_text,"fontsize",12);

figure;
plot([1:length(s_07)],imag(s_07),"b");
hold on;
plot([1:length(s_07_DMA_ideal)],imag(s_07_DMA_ideal),"r");
plot([1:length(theo_)],imag(theo_),"k--");
hold off;
legend("imbalance","ideal","theo");
title_text = "channel received compare (imag), practice 01  ";
title(title_text,"fontsize",12);
%}

s_10_normalized = real(s_10)/(sqrt(mean(real(s_10).^2)))+...
    imag(s_10)/(sqrt(mean(imag(s_10).^2))).*1i;

figure;
stem([1:length(s_10_normalized)],real(s_10_normalized),"bo");
hold on;
stem([1:length(s_1)],real(s_1),"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(real), practice 01  ";
title(title_text,"fontsize",12);

figure;
stem([1:length(s_10_normalized)],imag(s_10_normalized),"bo");
hold on;
stem([1:length(s_1)],imag(s_1),"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(imag), practice 01 ";
title(title_text,"fontsize",12);

%% pratice 02
% transmiter from LAB9
signal_length = 1e2;
s_1 = sign(randn(1,signal_length))+1i*sign(randn(1,signal_length));

H = [1 -imbalance_g*sin(imbalance_f) ; 0  imbalance_g*cos(imbalance_f)];

s_1_b = inv(H) * [real(s_1);imag(s_1)];
s_1_compensate = s_1_b(1,:)+s_1_b(2,:)*1i;
 
% upsampling factor is 4
s_01_up_1 = up_sample(M_1,s_1_compensate); % rate = symbol_rate*M

% SRRC
s_02_SRRC_1 = conv(s_01_up_1,SRRC_1_n);
s_02_SRRC_1 = s_02_SRRC_1([floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+1 :...
            floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning
        
s_03_up_2 = up_sample(M_2,s_02_SRRC_1);       

s_04_IIR = filter(Lab9_demo_IIR_2,s_03_up_2);

%s_04_IIR = filter(Lab9_demo_IIR_2,s_04_IIR);
%=============================================
% signal formation IQ imbalance lab10 p.3
xb_real = real(s_04_IIR);
xb_imag = imag(s_04_IIR);

fc = 1/4*DMA_rate; 

x_ideal = xb_real.*sqrt(2).*cos(2*pi*fc.*[1:length(xb_real)]/DMA_rate)-...
    xb_imag.*sqrt(2).*sin(2*pi*fc.*[1:length(xb_imag)]/DMA_rate);

x = xb_real.*sqrt(2).*cos(2*pi*fc.*[1:length(xb_real)]/DMA_rate)-...
    imbalance_g.*xb_imag.*sqrt(2).*sin(2*pi*fc.*[1:length(xb_imag)]/DMA_rate + imbalance_f);
%=============================================
carrier_rece = exp(1i.*2*pi.*fc.*[1:length(x)]/DMA_rate);

s_06_rece = x ./ carrier_rece;
s_07_DMA = filter(Lab9_demo_IIR_2,s_06_rece);

s_06_rece_ideal = x_ideal ./ carrier_rece;
s_07_DMA_ideal= filter(Lab9_demo_IIR_2,s_06_rece_ideal);

figure;
plot([1:length(s_07_DMA)],real(s_07_DMA),"b");
hold on;
plot([1:length(s_07_DMA_ideal)],real(s_07_DMA_ideal),"r--");
hold off;
legend("compensate","ideal");
title_text = "compensate received compare (real), practice 02  ";
title(title_text,"fontsize",12);

figure;
plot([1:length(s_07_DMA)],imag(s_07_DMA),"b");
hold on;
plot([1:length(s_07_DMA_ideal)],imag(s_07_DMA_ideal),"r--");
hold off;
legend("compensate","ideal");
title_text = "compensate received compare (imag), practice 02  ";
title(title_text,"fontsize",12);

delay = 8;
s_08_reDig_p1 = s_07_DMA([mod(delay,M_re_p1)+1:M_re_p1:length(s_07_DMA)]);

s_08_reDig_p1_normalized = real(s_08_reDig_p1)/(sqrt(mean(real(s_08_reDig_p1).^2)))+...
    imag(s_08_reDig_p1)/(sqrt(mean(imag(s_08_reDig_p1).^2))).*1i;

figure;
stem([1:length(s_08_reDig_p1_normalized)],real(s_08_reDig_p1_normalized),"bo");
hold on;
stem([1:length(s_1)],real(s_1),"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(real), practice 02  ";
title(title_text,"fontsize",12);

figure;
stem([1:length(s_08_reDig_p1_normalized)],imag(s_08_reDig_p1_normalized),"bo");
hold on;
stem([1:length(s_1)],imag(s_1),"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(imag), practice 02 ";
title(title_text,"fontsize",12);