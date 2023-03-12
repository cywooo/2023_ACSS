clear all; clc; close all;
%% series
num_of_series = 12000;
m = (sign(randn(1,num_of_series))+1)/2;

%% PAM

R_of_PAM = 3 ; 
N_of_PAM = 2^R_of_PAM ;

desicion = [0:N_of_PAM-1]*2-(N_of_PAM-1);
SNR_PAM_dB = 10;
SNR_PAM = 10^(SNR_PAM_dB/10);
PAM_mo_power = mean(desicion.^2);

noise_PAM_power = PAM_mo_power/SNR_PAM;

m_PAM_modu = reshape(m,[R_of_PAM,length(m)/R_of_PAM])';
m_PAM_modu = (bin2dec(num2str(m_PAM_modu))'*2)-(N_of_PAM-1);

%figure(1);
%stem([1:length(m_PAM_modu)],m_PAM_modu);

m_PAM_ch = m_PAM_modu + randn(1,length(m_PAM_modu))*sqrt(noise_PAM_power) ;

% de
for a = 1:(num_of_series/R_of_PAM)
    aa(a) = find(abs(desicion - m_PAM_ch(a)) == min(abs(desicion - m_PAM_ch(a))))-1;
end

Symbol_Error_PAM = sum(1.*(m_PAM_modu ~= aa))/num_of_series;

% ASII code of '0' is 48
aa = zeros(length(m)/R_of_PAM,R_of_PAM) + dec2bin(aa) -48;

m_demod_PAM = reshape(aa',[num_of_series,1])';

Bit_Error_PAM = sum(abs(m_demod_PAM-m))/num_of_series;


% certification
test_SNR_PAM = log10(mean(m_PAM_modu.^2)/mean((randn(1,length(m_PAM_modu))*sqrt(noise_PAM_power)).^2))*10;

%% QAM
N_of_QAM = 16;
layout = [2 2];

desicion_I = [0:2^layout(2)-1]*2-(2^layout(2)-1);
desicion_Q = [0:2^layout(1)-1]*2-(2^layout(1)-1);
graycode_I = [[0 0] ; [0 1] ; [1 1] ; [1 0]];
graycode_Q = [[0 0] ; [0 1] ; [1 1] ; [1 0]];

SNR_QAM_dB = 10;
SNR_QAM = 10^(SNR_QAM_dB/10);

%%%%%
QAM_mo_power = (18*4+10*8+4)/16;

noise_QAM_power = QAM_mo_power/SNR_QAM;

R_of_QAM = log2(N_of_QAM);
div_by_R = reshape(m,[R_of_QAM,length(m)/R_of_QAM])';

% 00->-3 01->-1 11-> 1 10-> 3
m_I_phase_bin = div_by_R(:,[1:layout(2)]);
m_Q_phase_bin = div_by_R(:,[layout(2)+1:R_of_QAM]);

%m_I_phase_Amp = (bin2dec(num2str(m_I_phase_bin))'*2)-(2^layout(2)-1) ;
%m_Q_phase_Amp = (bin2dec(num2str(m_Q_phase_bin))'*2)-(2^layout(1)-1) ;

% mo
for b = 1:(num_of_series/R_of_QAM)
    I_desi = graycode_I - m_I_phase_bin(b,:);
    Q_desi = graycode_Q - m_Q_phase_bin(b,:);
    I = abs(I_desi(:,1))+abs(I_desi(:,2));
    Q = abs(Q_desi(:,1))+abs(Q_desi(:,2));
    m_I_phase_Amp(b) = (find(I' == 0)-1)*2-(2^layout(2)-1); 
    m_Q_phase_Amp(b) = (find(Q' == 0)-1)*2-(2^layout(1)-1);
end

m_I_phase_ch = m_I_phase_Amp + randn(1,length(m_I_phase_Amp))*sqrt(noise_QAM_power) ;
m_Q_phase_ch = m_Q_phase_Amp + randn(1,length(m_I_phase_Amp))*sqrt(noise_QAM_power) ;

m_I_phase_demo = zeros(length(m)/R_of_QAM ,layout(2));
m_Q_phase_demo = zeros(length(m)/R_of_QAM ,layout(1));

% de
for b = 1:(num_of_series/R_of_QAM)
    bb_I(b) = find(abs(desicion_I - m_I_phase_ch(b)) == min(abs(desicion_I - m_I_phase_ch(b))))-1;
    m_I_phase_demo(b,:) = graycode_I(bb_I(b)+1,:);
    bb_Q(b) = find(abs(desicion_Q - m_Q_phase_ch(b)) == min(abs(desicion_Q - m_Q_phase_ch(b))))-1;
    m_Q_phase_demo(b,:) = graycode_Q(bb_Q(b)+1,:);
end

Symbol_Error_QAM = 1-sum(1.*(m_I_phase_Amp == bb_I*2-(2^layout(2)-1)).*(m_Q_phase_Amp == bb_Q*2-(2^layout(1)-1)))/num_of_series;

m_demo_QAM = [m_I_phase_demo m_Q_phase_demo];
m_demo_QAM = reshape(m_demo_QAM',[1,length(m)]);



figure(2);
scatter(m_I_phase_ch,m_Q_phase_ch);

Bit_Error_QAM = sum(abs(m_demo_QAM-m))/num_of_series;

% certification
test_SNR_QAM = log10(mean(m_I_phase_Amp.^2+m_Q_phase_Amp.^2)/mean((randn(1,length(m_PAM_modu))*sqrt(noise_PAM_power)).^2)*2)*10;

