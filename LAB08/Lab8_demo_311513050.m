clear all; clc; close all;

M_1 = 4;
M_2 = 4;

n_1 = [-3*M_1:3*M_1] + 1e-6;% Avoid Singularity
n_2 = [-3*M_2:3*M_2] + 1e-6;% Avoid Singularity

%% Pulse fomula SRRC
a = 0.5+1e-6 ;% Avoid Singularity

% n domain SRRC
SRRC_1_n = (4.*a./pi).*(cos((1+a).*pi.*n_1./M_1)+ M_1.*sin((1-a).*pi.*n_1./M_1)./(4.*a.*n_1))...
            ./(1-(4.*a.*n_1./M_1).^2);
        
SRRC_2_n = (4.*a./pi).*(cos((1+a).*pi.*n_2./M_2)+ M_2.*sin((1-a).*pi.*n_2./M_2)./(4.*a.*n_2))...
            ./(1-(4.*a.*n_2./M_2).^2);
        
%% practice 1
% 1.Conduct SRRC pulse shaping for a QPSK sequence (the upsampling factor is 64).
bpsk_length = 1e2;
s_BPSK = sign(randn(1,bpsk_length));


% upsampling factor is 4
s_01_up_1 = up_sample(M_1,s_BPSK); % rate = symbol_rate*M

% SRRC
s_02_SRRC_1 = conv(s_01_up_1,SRRC_1_n);
s_02_SRRC_1 = s_02_SRRC_1([floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+1 :...
            floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning

s_03_up_2 = up_sample(M_2,s_02_SRRC_1);       

s_04_SRRC_2 = conv(s_03_up_2,SRRC_2_n);
s_04_SRRC_2 = s_04_SRRC_2([floor((length(s_04_SRRC_2)-length(s_03_up_2))/2)+1 :...
            floor((length(s_04_SRRC_2)-length(s_03_up_2))/2)+length(s_03_up_2)]);% aligning
     
s_out = conv(s_04_SRRC_2,SRRC_2_n);
s_out = s_out([floor((length(s_out)-length(s_04_SRRC_2))/2)+1 :...
            floor((length(s_out)-length(s_04_SRRC_2))/2)+length(s_04_SRRC_2)]);% aligning

        
%% practice 2
s_05_add_noise = s_04_SRRC_2 + randn(1,length(s_04_SRRC_2)).*sqrt(0);

s_06_received_SRRC2 = conv(s_05_add_noise,SRRC_2_n);
s_06_received_SRRC2 = s_06_received_SRRC2([floor((length(s_06_received_SRRC2)-length(s_05_add_noise))/2)+1 :...
            floor((length(s_06_received_SRRC2)-length(s_05_add_noise))/2)+length(s_05_add_noise)]);% aligning

s_07_down_1 = down_sample(M_2,s_06_received_SRRC2);

s_08_received_SRRC1 = conv(s_07_down_1,SRRC_1_n);
s_08_received_SRRC1 = s_08_received_SRRC1([floor((length(s_08_received_SRRC1)-length(s_07_down_1))/2)+1 :...
            floor((length(s_08_received_SRRC1)-length(s_07_down_1))/2)+length(s_07_down_1)]);% aligning

a_hat_09 = down_sample(M_1,s_08_received_SRRC1);

a_hat_normalized = a_hat_09/(sqrt(mean(a_hat_09.^2)));

figure;
stem([1:length(a_hat_normalized)],real(a_hat_normalized));
hold on;
stem([1:length(s_BPSK)],real(s_BPSK));
hold off;
legend("a hat","s BPSK");

%% practice 3
% different since s_03_up_2 

s_04_IIR = filter(Lab8_demo_IIR,s_03_up_2);

s_05_IIR_add_noise = s_04_IIR + randn(1,length(s_04_IIR)).*sqrt(0);

s_06_IIR_received = filter(Lab8_demo_IIR,s_05_IIR_add_noise);

s_07_IIR_down_1 = down_sample(M_2,s_06_IIR_received);

s_08_IIR_SRRC1 = conv(s_07_IIR_down_1,SRRC_1_n);
s_08_IIR_SRRC1 = s_08_IIR_SRRC1([floor((length(s_08_IIR_SRRC1)-length(s_07_IIR_down_1))/2)+1 :...
            floor((length(s_08_IIR_SRRC1)-length(s_07_IIR_down_1))/2)+length(s_07_IIR_down_1)]);% aligning

a_hat_IIR_09 = down_sample(M_1,s_08_IIR_SRRC1);

a_hat_IIR_normalized = a_hat_IIR_09/(sqrt(mean(a_hat_IIR_09.^2)));

figure;
stem([1:length(a_hat_IIR_normalized)],real(a_hat_IIR_normalized));
hold on;
stem([1:length(s_BPSK)],real(s_BPSK));
hold off;
legend("a hat IIR","s BPSK");
title_text = "p3";
title(title_text,"fontsize",12);

