clear all; clc; close all;
%% series
num_of_series = 1.2e4;
m = (sign(randn(1,num_of_series))+1)/2;

%% PAM

PAM_R = 3 ; 
PAM_N = 2^PAM_R ;

PAM_desicion = [0:PAM_N-1]*2-(PAM_N-1);

%%% SNR adjust
PAM_SNR_dB = 10;
%%% SNR adjust

PAM_SNR = 10^(PAM_SNR_dB/10);
PAM_mo_power = mean(PAM_desicion.^2);

PAM_noise_power = PAM_mo_power/PAM_SNR;

PAM_m_modu = reshape(m,[PAM_R,length(m)/PAM_R])';
PAM_m_modu = (bin2dec(num2str(PAM_m_modu))'*2)-(PAM_N-1);

PAM_m_ch = PAM_m_modu + randn(1,length(PAM_m_modu))*sqrt(PAM_noise_power) ;

% de
for a = 1:(num_of_series/PAM_R)
    aa(a) = find(abs(PAM_desicion - PAM_m_ch(a)) == min(abs(PAM_desicion - PAM_m_ch(a))))-1;
end

PAM_Error_Symbol = sum(1.*(PAM_m_modu ~= aa))/num_of_series;

% ASII code of '0' is 48
aa = zeros(length(m)/PAM_R,PAM_R) + dec2bin(aa) -48;

PAM_m_demod = reshape(aa',[num_of_series,1])';

PAM_Error_Bit = sum(abs(PAM_m_demod-m))/num_of_series;

figure(1);
plot([1:length(PAM_m_modu)],PAM_m_modu);

% certification
PAM_SNR_test = log10(mean(PAM_m_modu.^2)/mean((randn(1,length(PAM_m_modu))*sqrt(PAM_noise_power)).^2))*10;

%% QAM
QAM_N = 16;
QAM_R = log2(QAM_N);
layout = [2 2];

QAM_desicion_I = [0:2^layout(2)-1]*2-(2^layout(2)-1);
QAM_desicion_Q = [0:2^layout(1)-1]*2-(2^layout(1)-1);
QAM_graycode_I = [[0 0] ; [0 1] ; [1 1] ; [1 0]];
QAM_graycode_Q = [[0 0] ; [0 1] ; [1 1] ; [1 0]];

%%% SNR adjust
QAM_SNR_dB = 10;
%%% SNR adjust

QAM_SNR = 10^(QAM_SNR_dB/10);

%%%%%
QAM_mo_power = 0;
for k = 1:layout(1)^2
    for q = 1:layout(2)^2
        QAM_mo_power = QAM_mo_power + QAM_desicion_I(k)^2 + QAM_desicion_Q(q)^2 ;
    end
end
QAM_mo_power = QAM_mo_power/QAM_N;

QAM_noise_power = QAM_mo_power/QAM_SNR;

QAM_div_by_R = reshape(m,[QAM_R,length(m)/QAM_R])';

% 00->-3 01->-1 11-> 1 10-> 3
QAM_m_I_phase_bin = QAM_div_by_R(:,[1:layout(2)]);
QAM_m_Q_phase_bin = QAM_div_by_R(:,[layout(2)+1:QAM_R]);

% mo
for b = 1:(num_of_series/QAM_R)
    I_desi = QAM_graycode_I - QAM_m_I_phase_bin(b,:);
    Q_desi = QAM_graycode_Q - QAM_m_Q_phase_bin(b,:);
    I = abs(I_desi(:,1))+abs(I_desi(:,2));
    Q = abs(Q_desi(:,1))+abs(Q_desi(:,2));
    m_I_phase_Amp(b) = (find(I' == 0)-1)*2-(2^layout(2)-1); 
    m_Q_phase_Amp(b) = (find(Q' == 0)-1)*2-(2^layout(1)-1);
end

QAM_m_I_phase_ch = m_I_phase_Amp + randn(1,length(m_I_phase_Amp))*sqrt(QAM_noise_power/2) ;
QAM_m_Q_phase_ch = m_Q_phase_Amp + randn(1,length(m_Q_phase_Amp))*sqrt(QAM_noise_power/2) ;

QAM_m_I_phase_demo = zeros(length(m)/QAM_R ,layout(2));
QAM_m_Q_phase_demo = zeros(length(m)/QAM_R ,layout(1));

% de
for b = 1:(num_of_series/QAM_R)
    bb_I(b) = find(abs(QAM_desicion_I - QAM_m_I_phase_ch(b)) == min(abs(QAM_desicion_I - QAM_m_I_phase_ch(b))));
    QAM_m_I_phase_demo(b,:) = QAM_graycode_I(bb_I(b),:);
    bb_Q(b) = find(abs(QAM_desicion_Q - QAM_m_Q_phase_ch(b)) == min(abs(QAM_desicion_Q - QAM_m_Q_phase_ch(b))));
    QAM_m_Q_phase_demo(b,:) = QAM_graycode_Q(bb_Q(b),:);
end

QAM_count = sum(abs([QAM_m_I_phase_demo QAM_m_Q_phase_demo]-[QAM_m_I_phase_bin QAM_m_Q_phase_bin])');
QAM_Error_Symbol = 1-sum(1.*(QAM_count == 0))/length(QAM_count);

QAM_m_demo = [QAM_m_I_phase_demo QAM_m_Q_phase_demo];
QAM_m_demo = reshape(QAM_m_demo',[1,length(m)]);



figure(2);
scatter(QAM_m_I_phase_ch,QAM_m_Q_phase_ch);

QAM_Error_Bit = sum(abs(QAM_m_demo-m))/ ;

% certification
QAM_test_SNR = log10(mean(m_I_phase_Amp.^2+m_Q_phase_Amp.^2)/QAM_noise_power^2)*10;

QAM_p = 3/4*erfc(sqrt(QAM_SNR/10));
QAM_Error_Symbol_theory = 1- (1-QAM_p)^2;
