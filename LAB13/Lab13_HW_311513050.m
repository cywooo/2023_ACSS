clear all; clc; close all;

symbol_rate = 1e6;
DAC_rate = 16e6;
DMA_rate = 32e6;
analog_freq = 8e6;

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

SNR = [-20:5:10]; % inf;%    
guard_bit = 20;
batch = 500;
%% 
for q = [1:length(SNR)]
    for qq = [1:batch]
        signal_length = 1e3 + guard_bit; % signal length adjust
        % signal generate------
        s_2 = [zeros(1,guard_bit/2) sign(randn(1,signal_length-guard_bit)) zeros(1,guard_bit/2)]; %
        % up sample by 16-------
        s_01_up_1 = up_sample(M_1,s_2);
        % SRRC filter ---------
        s_02_SRRC_1 = conv(s_01_up_1,SRRC_1_n);
        s_02_SRRC_1 = s_02_SRRC_1([floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+1 :...
                    floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning           
        %-----------------------
        % up sample by 2-------
        s_03_up_2 = up_sample(M_2,s_02_SRRC_1);   
        % DMA filter------------
        s_04_IIR = conv(s_03_up_2,SRRC_2_n);
        s_04_IIR = s_04_IIR([floor((length(s_04_IIR)-length(s_03_up_2))/2)+1 :...
                    floor((length(s_04_IIR)-length(s_03_up_2))/2)+length(s_03_up_2)]);% aligning
        %-----------------------
        % upconvert----------
        fc = 1/4*DMA_rate;
        x =  s_04_IIR.*exp(1i*2*pi*fc.*[1:length(s_04_IIR)]./DMA_rate);
        x = real(x);
        % Channel===============       
        h_2_to_2 = [ randn(1,1)*exp(-j*2*pi*rand(1,1)) randn(1,1)*exp(-j*2*pi*rand(1,1)) ;...
                     randn(1,1)*exp(-j*2*pi*rand(1,1)) randn(1,1)*exp(-j*2*pi*rand(1,1)) ];         
        
        num_of_antana = 2;       
        y_2_to_2 = h_2_to_2* [ x ;...
                               x ];
        
        [y_x,y_y] = size(y_2_to_2);
        
        noise_power = 1/(10^(SNR(q)/10));
        
        y_channel = y_2_to_2 + randn(y_x,y_y)*sqrt(noise_power);
        
        W_ZF = inv(h_2_to_2);
        y_ZF = W_ZF * y_channel;
        
        W_MMSE = inv( (h_2_to_2') * h_2_to_2 + (10^(SNR(q)/10))^-1 *eye(num_of_antana)) * (h_2_to_2');
        y_MMSE = W_MMSE * y_channel;
        
        % downconvert---------
        s_06_ZF_rece = y_ZF .*exp(-1*1i*2*pi*fc.*[1:length(y_ZF)]./DMA_rate);
        s_06_MMSE_rece = y_MMSE .*exp(-1*1i*2*pi*fc.*[1:length(y_MMSE)]./DMA_rate);
        % received DMA------------------        
        s_07_ZF_DMA_split = [];
        for k=[1:num_of_antana]        
            s_07_ZF_DMA_split(k,:) = conv(s_06_ZF_rece(k,:),SRRC_2_n);
        end
        s_07_ZF_DMA_split = s_07_ZF_DMA_split(:,[floor((length(s_07_ZF_DMA_split)-length(s_06_ZF_rece))/2)+1 :...
                    floor((length(s_07_ZF_DMA_split)-length(s_06_ZF_rece))/2)+length(s_06_ZF_rece)]);% aligning
        
        s_07_MMSE_DMA_split = [];
        for k=[1:num_of_antana]        
            s_07_MMSE_DMA_split(k,:) = conv(s_06_MMSE_rece(k,:),SRRC_2_n);
        end
        s_07_MMSE_DMA_split = s_07_MMSE_DMA_split(:,[floor((length(s_07_MMSE_DMA_split)-length(s_06_MMSE_rece))/2)+1 :...
                    floor((length(s_07_MMSE_DMA_split)-length(s_06_MMSE_rece))/2)+length(s_06_MMSE_rece)]);% aligning
        % down sample by 2--------------
        s_08_ZF_reDig_split = s_07_ZF_DMA_split(:,[1:M_2:length(s_07_ZF_DMA_split)]);
        s_09_ZF_toDigFilter_split = s_08_ZF_reDig_split;
        
        s_08_MMSE_reDig_split = s_07_MMSE_DMA_split(:,[1:M_2:length(s_07_MMSE_DMA_split)]);
        s_09_MMSE_toDigFilter_split = s_08_MMSE_reDig_split;
        % SRRC filter ------------------
        s_10_ZF_DigFilter_split=[];
        for k=[1:num_of_antana]      
            s_10_ZF_DigFilter_split(k,:) = conv(s_09_ZF_toDigFilter_split(k,:),SRRC_1_n);
        end
        s_10_ZF_DigFilter_split = s_10_ZF_DigFilter_split(:,[floor((length(s_10_ZF_DigFilter_split)-length(s_09_ZF_toDigFilter_split))/2)+1 :...
                    floor((length(s_10_ZF_DigFilter_split)-length(s_09_ZF_toDigFilter_split))/2)+length(s_09_ZF_toDigFilter_split)]);% aligning      
        
        s_10_MMSE_DigFilter_split=[];
        for k=[1:num_of_antana]      
            s_10_MMSE_DigFilter_split(k,:) = conv(s_09_MMSE_toDigFilter_split(k,:),SRRC_1_n);
        end
        s_10_MMSE_DigFilter_split = s_10_MMSE_DigFilter_split(:,[floor((length(s_10_MMSE_DigFilter_split)-length(s_09_MMSE_toDigFilter_split))/2)+1 :...
                    floor((length(s_10_MMSE_DigFilter_split)-length(s_09_MMSE_toDigFilter_split))/2)+length(s_09_MMSE_toDigFilter_split)]);% aligning      
        
        %-------------------------------
        % down sample by 16-------------
        delay = 1;
        s_11_ZF_reDig_split = s_10_ZF_DigFilter_split(:,[delay:M_1:length(s_10_ZF_DigFilter_split)]);
        s_11_ZF_reDig = sum(s_11_ZF_reDig_split)/32;
        
        %s_11_ZF_reDig_normalized = s_11_ZF_reDig/sqrt(mean(abs(s_11_ZF_reDig([guard_bit/2+1:end-guard_bit/2])).^2));
        s_11_ZF_reDig_normalized = s_11_ZF_reDig([guard_bit/2+1:end-guard_bit/2]);
        
        s_11_MMSE_reDig_split = s_10_MMSE_DigFilter_split(:,[delay:M_1:length(s_10_MMSE_DigFilter_split)]);
        s_11_MMSE_reDig = sum(s_11_MMSE_reDig_split);
        
        %s_11_MMSE_reDig_normalized = s_11_MMSE_reDig/sqrt(mean(abs(s_11_MMSE_reDig([6:end-5])).^2));
        s_11_MMSE_reDig_normalized = s_11_MMSE_reDig([guard_bit/2+1:end-guard_bit/2])/30;
        
        BER_ZF_matrix(q,qq) = mean(abs(sign(real(s_11_ZF_reDig_normalized))-s_2([guard_bit/2+1:end-guard_bit/2]))/2);
        
        BER_MMSE_matrix(q,qq) = mean(abs(sign(real(s_11_MMSE_reDig_normalized))-s_2([guard_bit/2+1:end-guard_bit/2]))/2);

        OSNR_ZF = 1/mean(abs(s_11_ZF_reDig_normalized-s_2([guard_bit/2+1:end-guard_bit/2])).^2);
        OSNR_ZF_batch_matrix(q,qq) = OSNR_ZF;
        
        OSNR_MMSE = 1/mean(abs(s_11_MMSE_reDig_normalized-s_2([guard_bit/2+1:end-guard_bit/2])).^2);
        OSNR_MMSE_batch_matrix(q,qq) = OSNR_MMSE;
    end
end

BER_ZF = mean(BER_ZF_matrix')';

BER_MMSE = mean(BER_MMSE_matrix')';

OSNR_ZF_batch = mean(OSNR_ZF_batch_matrix')';
OSNR_ZF_dB   = 10*log10(OSNR_ZF_batch);

OSNR_MMSE_batch = mean(OSNR_MMSE_batch_matrix')';
OSNR_MMSE_dB = 10*log10(OSNR_MMSE_batch);

figure(1);
semilogy(SNR,BER_ZF,"b");
hold on;
semilogy(SNR,BER_MMSE,"r");
hold off;
legend("ZF","MMSE");
xlabel("SNR");
ylabel("BER");

figure(2);
plot(SNR,OSNR_ZF_dB,"b");
hold on;
plot(SNR,OSNR_MMSE_dB,"r");
hold off;
legend("ZF","MMSE");
xlabel("SNR");
ylabel("SNR_o");