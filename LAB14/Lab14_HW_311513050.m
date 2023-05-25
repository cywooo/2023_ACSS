clear all; clc; close all;

% Use Lab09 structure

symbol_rate = 1e6;
DAC_rate = 16e6;
DMA_rate = 64e6;
carrier_freq_analog = 16e6;
IF_freq = 4e6;

M_1 = DAC_rate/symbol_rate;
M_2 = DMA_rate/DAC_rate;


n_1 = [-2*M_1:2*M_1] + 1e-6;% Avoid Singularity
n_2 = [-2*M_2:2*M_2] + 1e-6;% Avoid Singularity
a = 0.8+1e-6 ;% Avoid Singularity

% n domain SRRC     
A = cos((1+a).*pi.*n_1./M_1);
B = M_1.*sin((1-a).*pi.*n_1./M_1)./(4.*a.*n_1);
C = 1-(4.*a.*n_1./M_1).^2;
SRRC_1_n = (4.*a./pi).*(A+B)./C;
% max(SRRC_1_n) = 1.2186 & min(SRRC_1_n) = -0.1104
nob_SRRC_1 = 5;
DR_SRRC_1 = 2^1; %+-
code_map_SRRC_1 = -DR_SRRC_1:((DR_SRRC_1*2)/(2^nob_SRRC_1)):(DR_SRRC_1-(DR_SRRC_1*2)/(2^nob_SRRC_1));
SRRC_1_F = F_point_decision(SRRC_1_n,code_map_SRRC_1);

A = cos((1+a).*pi.*n_2./M_2);
B = M_2.*sin((1-a).*pi.*n_2./M_2)./(4.*a.*n_2);
C = 1-(4.*a.*n_2./M_2).^2;
SRRC_2_n = (4.*a./pi).*(A+B)./C;
% max(SRRC_2_n) = 1.2186 & min(SRRC_2_n) = -0.1094
nob_SRRC_2 = 5;
DR_SRRC_2 = 2^1; %+-
code_map_SRRC_1 = -DR_SRRC_2:((DR_SRRC_2*2)/(2^nob_SRRC_2)):(DR_SRRC_2-(DR_SRRC_2*2)/(2^nob_SRRC_2));
SRRC_2_F = F_point_decision(SRRC_2_n,code_map_SRRC_2);
        
bpsk_length = 40 + 10;
s_BPSK =[zeros(1,5) sign(randn(1,bpsk_length-10)) zeros(1,5)];
%% no Float
% upsampling factor is 16
s_01_up_1 = up_sample(M_1,s_BPSK); % rate = symbol_rate*M

% SRRC
s_02_SRRC_1 = conv(s_01_up_1,SRRC_1_n);
s_02_SRRC_1 = s_02_SRRC_1([floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+1 :...
            floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning
        
s_03_up_2 = up_sample(M_2,s_02_SRRC_1);       

s_04_SRRC_2 = conv(s_03_up_2,SRRC_2_n);
s_04_SRRC_2 = s_04_SRRC_2([floor((length(s_04_SRRC_2)-length(s_03_up_2))/2)+1 :...
            floor((length(s_04_SRRC_2)-length(s_03_up_2))/2)+length(s_03_up_2)]);% aligning

w_c = carrier_freq_analog/symbol_rate/M_1/M_2*2*pi;
carrier_tran = exp(i*w_c.*[1:length(s_04_SRRC_2)]);
carrier_rece = exp(-i*w_c.*[1:length(s_04_SRRC_2)]);

s_05_channel = real(s_04_SRRC_2 .* carrier_tran);
s_05_channel_noise = s_05_channel + randn(1,length(s_05_channel))*sqrt(0);

s_06_IF_rece = s_05_channel_noise .* carrier_rece;

s_07_IF_SRRC_2 = conv(s_06_IF_rece,SRRC_2_n);
s_07_IF_SRRC_2 = s_07_IF_SRRC_2([floor((length(s_07_IF_SRRC_2)-length(s_06_IF_rece))/2)+1 :...
            floor((length(s_07_IF_SRRC_2)-length(s_06_IF_rece))/2)+length(s_06_IF_rece)]);% aligning

delay = 1;
s_08_IF_reDig = s_07_IF_SRRC_2([mod(delay,M_2)+1:M_2:length(s_07_IF_SRRC_2)]);
 
s_09_IF_toDigFilter = s_08_IF_reDig;
s_10_IF_DigFilter = conv(s_09_IF_toDigFilter,SRRC_1_n);
s_10_IF_DigFilter = s_10_IF_DigFilter([floor((length(s_10_IF_DigFilter)-length(s_09_IF_toDigFilter))/2)+1 :...
            floor((length(s_10_IF_DigFilter)-length(s_09_IF_toDigFilter))/2)+length(s_09_IF_toDigFilter)]);% aligning

delay = 1;
s_11_IF_reDig = s_10_IF_DigFilter([mod(delay,M_1)+1:M_1:length(s_10_IF_DigFilter)]);
s_11_IF_reDig = s_11_IF_reDig([6:end-5]);

s_11_IF_reDig_normalized = s_11_IF_reDig/(sqrt(mean(real(s_11_IF_reDig).^2)));

SQNR = 10*log10(norm(s_BPSK([6:end-5]))/ norm(abs(s_11_IF_reDig_normalized-s_BPSK([6:end-5]))));
%% Floating
% upsampling factor is 16
sf_01_up_1 = up_sample(M_1,s_BPSK); % rate = symbol_rate*M

% SRRC
sf_02_SRRC_1 = conv(sf_01_up_1,SRRC_1_F);
sf_02_SRRC_1 = sf_02_SRRC_1([floor((length(sf_02_SRRC_1)-length(sf_01_up_1))/2)+1 :...
            floor((length(sf_02_SRRC_1)-length(sf_01_up_1))/2)+length(sf_01_up_1)]);% aligning
% max(sf_02_SRRC_1) = 1.5 & min(sf_02_SRRC_1) = -1.5
nob = 5;
DR = 2^1; %+-
code_map = -DR:((DR*2)/(2^nob)):(DR-(DR*2)/(2^nob));
SRRC_1_F = F_point_decision(sf_02_SRRC_1,code_map);
        
s_03_up_2 = up_sample(M_2,s_02_SRRC_1);       

s_04_SRRC_2 = conv(s_03_up_2,SRRC_2_n);
s_04_SRRC_2 = s_04_SRRC_2([floor((length(s_04_SRRC_2)-length(s_03_up_2))/2)+1 :...
            floor((length(s_04_SRRC_2)-length(s_03_up_2))/2)+length(s_03_up_2)]);% aligning

w_c = carrier_freq_analog/symbol_rate/M_1/M_2*2*pi;
carrier_tran = exp(i*w_c.*[1:length(s_04_SRRC_2)]);
carrier_rece = exp(-i*w_c.*[1:length(s_04_SRRC_2)]);

s_05_channel = real(s_04_SRRC_2 .* carrier_tran);
s_05_channel_noise = s_05_channel + randn(1,length(s_05_channel))*sqrt(0);

s_06_IF_rece = s_05_channel_noise .* carrier_rece;

s_07_IF_SRRC_2 = conv(s_06_IF_rece,SRRC_2_n);
s_07_IF_SRRC_2 = s_07_IF_SRRC_2([floor((length(s_07_IF_SRRC_2)-length(s_06_IF_rece))/2)+1 :...
            floor((length(s_07_IF_SRRC_2)-length(s_06_IF_rece))/2)+length(s_06_IF_rece)]);% aligning

delay = 1;
s_08_IF_reDig = s_07_IF_SRRC_2([mod(delay,M_2)+1:M_2:length(s_07_IF_SRRC_2)]);
 
s_09_IF_toDigFilter = s_08_IF_reDig;
s_10_IF_DigFilter = conv(s_09_IF_toDigFilter,SRRC_1_n);
s_10_IF_DigFilter = s_10_IF_DigFilter([floor((length(s_10_IF_DigFilter)-length(s_09_IF_toDigFilter))/2)+1 :...
            floor((length(s_10_IF_DigFilter)-length(s_09_IF_toDigFilter))/2)+length(s_09_IF_toDigFilter)]);% aligning

delay = 1;
s_11_IF_reDig = s_10_IF_DigFilter([mod(delay,M_1)+1:M_1:length(s_10_IF_DigFilter)]);
s_11_IF_reDig = s_11_IF_reDig([6:end-5]);

s_11_IF_reDig_normalized = s_11_IF_reDig/(sqrt(mean(real(s_11_IF_reDig).^2)));

SQNR = 10*log10(norm(s_BPSK([6:end-5]))/ norm(abs(s_11_IF_reDig_normalized-s_BPSK([6:end-5]))));
%% Figure
figure;
stem([1:length(s_11_IF_reDig_normalized)],real(s_11_IF_reDig_normalized),"rO");
hold on;
stem([1:length(s_BPSK)-10],real(s_BPSK([6:end-5])),"b+");
hold off;
legend("IF Receive","s BPSK, Transmit");
title_text = "Transmit and Receive, SQNR = "+num2str(SQNR);
title(title_text,"fontsize",12);