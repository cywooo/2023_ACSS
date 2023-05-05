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

%% HW
signal_length = 50;
% signal generate------
s_real = [zeros(1,10) 1 zeros(1,signal_length-10-1)]; %sign(randn(1,signal_length)); %
s_imag = sign(randn(1,signal_length)); %[-1*ones(1,signal_length/2) ones(1,signal_length/2)]; %
s_1 = s_real + s_imag.*1i;
% up sample by 16-------
s_01_up_1 = up_sample(M_1,s_1);
% SRRC filter ---------
s_02_SRRC_1 = conv(s_01_up_1,SRRC_1_n);
s_02_SRRC_1 = s_02_SRRC_1([floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+1 :...
            floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning
        
%s_02_SRRC_1_real = filter(Lab10_HW_IIR_16,s_01_up_1_real);         
%s_02_SRRC_1_imag = filter(Lab10_HW_IIR_16,s_01_up_1_imag);                 
%-----------------------
% up sample by 2-------
s_03_up_2 = up_sample(M_2,s_02_SRRC_1);   
% DMA filter------------
%s_04_IIR_real = filter(Lab10_demo_IIR,s_03_up_2_real);
%s_04_IIR_imag = filter(Lab10_demo_IIR,s_03_up_2_imag);
s_04_IIR = conv(s_03_up_2,SRRC_2_n);
s_04_IIR = s_04_IIR([floor((length(s_04_IIR)-length(s_03_up_2))/2)+1 :...
            floor((length(s_04_IIR)-length(s_03_up_2))/2)+length(s_03_up_2)]);% aligning
%-----------------------
% IF upconvert----------
fc = 1/4*DMA_rate;
x =  s_04_IIR.*exp(1i*2*pi*fc.*[1:length(s_04_IIR)]./DMA_rate);
x = real(x);
% Channel---------------
M = M_1*M_2;
apha_b = [1 -.5];
tou = [0 1];
num_of_multipath = length(apha_b);
x_ = zeros(num_of_multipath ,length(x));
for k=[1:num_of_multipath]
    x_(k,:) = apha_b(k).*[zeros(1,floor(tou(k)*M)) x(1:(length(x)-floor(tou(k)*M)))];
end


y = sum(x_);

% IF downconvert to IF freq.----
wc = 2*pi*fc/DMA_rate;
wif = 2*pi*IF_freq/DMA_rate;

s_06_rece = y .*exp(-1i*2*pi*fc.*[1:length(s_04_IIR)]./DMA_rate);
s_06_rece_split = x_ .*exp(-1i*2*pi*fc.*[1:length(s_04_IIR)]./DMA_rate);
% received DMA------------------
%s_07_DMA = filter(Lab10_demo_IIR,s_06_rece); 

s_07_DMA = conv(s_06_rece,SRRC_2_n);
s_07_DMA = s_07_DMA([floor((length(s_07_DMA)-length(s_06_rece))/2)+1 :...
            floor((length(s_07_DMA)-length(s_06_rece))/2)+length(s_06_rece)]);% aligning
        
for k=[1:num_of_multipath]        
    s_07_DMA_split(k,:) = conv(s_06_rece_split(k,:),SRRC_2_n);
end
s_07_DMA_split = s_07_DMA_split(:,[floor((length(s_07_DMA_split)-length(s_06_rece_split))/2)+1 :...
            floor((length(s_07_DMA_split)-length(s_06_rece_split))/2)+length(s_06_rece_split)]);% aligning
% down sample by 2--------------
s_08_reDig = s_07_DMA([1:M_2:length(s_07_DMA)]);
s_09_IF_toDigFilter = s_08_reDig;

s_08_reDig_split = s_07_DMA_split(:,[1:M_2:length(s_07_DMA_split)]);
s_09_IF_toDigFilter_split = s_08_reDig_split;
% SRRC filter ------------------
s_10_IF_DigFilter = conv(s_09_IF_toDigFilter,SRRC_1_n);
s_10_IF_DigFilter = s_10_IF_DigFilter([floor((length(s_10_IF_DigFilter)-length(s_09_IF_toDigFilter))/2)+1 :...
            floor((length(s_10_IF_DigFilter)-length(s_09_IF_toDigFilter))/2)+length(s_09_IF_toDigFilter)]);% aligning

for k=[1:num_of_multipath]      
    s_10_IF_DigFilter_split(k,:) = conv(s_09_IF_toDigFilter_split(k,:),SRRC_1_n);
end
s_10_IF_DigFilter_split = s_10_IF_DigFilter_split(:,[floor((length(s_10_IF_DigFilter_split)-length(s_09_IF_toDigFilter_split))/2)+1 :...
            floor((length(s_10_IF_DigFilter_split)-length(s_09_IF_toDigFilter_split))/2)+length(s_09_IF_toDigFilter_split)]);% aligning
%s_10_IF_DigFilter = filter(Lab10_HW_IIR_16,s_09_IF_toDigFilter);         
%-------------------------------
% down sample by 16-------------
delay = 1;
s_11_IF_reDig = s_10_IF_DigFilter([delay:M_1:length(s_10_IF_DigFilter)]);
s_11_IF_reDig_split = s_10_IF_DigFilter_split(:,[delay:M_1:length(s_10_IF_DigFilter_split)]);
%-------------------------------
% normalized-------------
adj = 32/2;
s_11_IF_reDig_normalized = real(s_11_IF_reDig)/adj+...
    imag(s_11_IF_reDig)/adj*1i;
%(sqrt(mean(real(s_11_IF_reDig_split(k,:)).^2)))
for k=[1:num_of_multipath]     
    s_11_IF_reDig_split_normalized(k,:) = real(s_11_IF_reDig_split(k,:))/adj+...
        imag(s_11_IF_reDig_split(k,:))/adj*1i;
end
% equalizer--------------
equizer = eq(apha_b,tou,10);
equizer_ = equizer/sqrt(mean(equizer.^2));
equal_bit = conv(s_11_IF_reDig_normalized,equizer_);




figure(1);
subplot(2,1,1);
stem([1:signal_length],real(s_11_IF_reDig_normalized),"bo");
hold on;
stem([1:signal_length],s_real,"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(real), HW";
title(title_text,"fontsize",12);

subplot(2,1,2);
stem([1:signal_length],imag(s_11_IF_reDig_normalized),"bo");
hold on;
stem([1:signal_length],s_imag,"r+");
hold off;
legend("re","tr");
title_text = "Transmit and Receive(imag), HW";
title(title_text,"fontsize",12);

figure(2);
subplot(2,1,1);
stem([1:signal_length],real(s_11_IF_reDig_normalized),"bo");
hold on;
stem([1:signal_length],real(s_11_IF_reDig_split_normalized(1,:)),"r+");
stem([1:signal_length],real(s_11_IF_reDig_split_normalized(2,:)),"rX");
hold off;
legend("re","path 1","path 2");
title_text = "Receive and Path ingrident(real), HW";
title(title_text,"fontsize",12);

subplot(2,1,2);
stem([1:signal_length],imag(s_11_IF_reDig_normalized),"bo");
hold on;
stem([1:signal_length],imag(s_11_IF_reDig_split_normalized(1,:)),"r+");
stem([1:signal_length],imag(s_11_IF_reDig_split_normalized(2,:)),"rX");
hold off;
legend("re","path 1","path 2");
title_text = "Receive and Path ingrident(imag), HW";
title(title_text,"fontsize",12);

figure(4);
stem([1:length(equizer)],equizer,"bo");
title_text = "Equalizer";
title(title_text,"fontsize",12);

figure(3);
subplot(2,1,1);
stem([1:length(equal_bit)],real(equal_bit),"bo");
hold on;
stem([1:signal_length],s_real,"r+");
hold off;
legend("re_eq","tr");
title_text = "Transmit and Receive EQ (real), HW";
title(title_text,"fontsize",12);
subplot(2,1,2);
stem([1:length(equal_bit)],imag(equal_bit),"bo");
hold on;
stem([1:signal_length],s_imag,"r+");
hold off;
legend("re_eq","tr");
title_text = "Transmit and Receive EQ (imag), HW";
title(title_text,"fontsize",12);

function equalizer = eq(apha_b,tou,len)
    syms z n;
    %for 2 tap
    F = 0;
    for k = [1:length(apha_b)]      
        F = F + apha_b(k)*z^(-tou(k))
    end
    F^-1
    %equalizer = double(subs(sol,n,[0:len-1]))
    %{%}
    if(-apha_b(2)/apha_b(1))<=1
        sol = iztrans(F^-1)
        equalizer = double(subs(sol,n,[0:len-1]))
    end
    if(-apha_b(2)/apha_b(1))>1
        sol = -1*iztrans(F^-1)
        equalizer = double(subs(sol,n,-[0:len-1]-1))
    end
    
end
