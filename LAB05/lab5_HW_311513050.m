clear all; clc; close all;
%% series
num_of_series = 1.2e6;
m = (sign(randn(1,num_of_series))+1)/2;

%% QAM
QAM_N = 16;
QAM_R = log2(QAM_N);
layout = [2 2];

QAM_desicion_I = [0:2^layout(2)-1]*2-(2^layout(2)-1);
QAM_desicion_Q = [0:2^layout(1)-1]*2-(2^layout(1)-1);
QAM_graycode_I = [[0 0] ; [0 1] ; [1 1] ; [1 0]];
QAM_graycode_Q = [[0 0] ; [0 1] ; [1 1] ; [1 0]];


%%% SNR adjust
%QAM_SNR_dB = 10;
%%% SNR adjust
QAM_SNR_dB = [1:20];
show_fig = [5 10 15 20];
QAM_Error_Symbol_theory = zeros(1,length(QAM_SNR_dB));
QAM_Error_Symbol = zeros(1,length(QAM_SNR_dB));
a=1;
for g = 1:length(QAM_SNR_dB);
    QAM_SNR = 10^(QAM_SNR_dB(g)/10);
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
    for b = 1:(num_of_series/QAM_R);
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
    for b = 1:(num_of_series/QAM_R);
        bb_I(b) = find(abs(QAM_desicion_I - QAM_m_I_phase_ch(b)) == min(abs(QAM_desicion_I - QAM_m_I_phase_ch(b))));
        QAM_m_I_phase_demo(b,:) = QAM_graycode_I(bb_I(b),:);
        bb_Q(b) = find(abs(QAM_desicion_Q - QAM_m_Q_phase_ch(b)) == min(abs(QAM_desicion_Q - QAM_m_Q_phase_ch(b))));
        QAM_m_Q_phase_demo(b,:) = QAM_graycode_Q(bb_Q(b),:);
    end
    
    QAM_count = sum(abs([QAM_m_I_phase_demo QAM_m_Q_phase_demo]-[QAM_m_I_phase_bin QAM_m_Q_phase_bin])');
    QAM_Error_Symbol(g) = 1-sum(1.*(QAM_count == 0))/length(QAM_count);
    
    QAM_m_demo = [QAM_m_I_phase_demo QAM_m_Q_phase_demo];
    QAM_m_demo = reshape(QAM_m_demo',[1,length(m)]);
    
    if max(QAM_SNR_dB(g) == show_fig);
        figure(1);
        subplot(2,2,a);
        scatter(QAM_m_I_phase_ch,QAM_m_Q_phase_ch,2);
        axis([-8 8 -8 8]);
        title("SNR = ",QAM_SNR_dB(g));
        a = a+1;
    end
    
    QAM_Error_Bit = sum(abs(QAM_m_demo-m))/num_of_series;
    
    % certification
    QAM_test_SNR = log10(mean(m_I_phase_Amp.^2+m_Q_phase_Amp.^2)/QAM_noise_power^2)*10;
    
    QAM_p = 3/4*erfc(sqrt(QAM_SNR/10));
    QAM_Error_Symbol_theory(g) = 1- (1-QAM_p)^2;
end

figure;
semilogy([1:length(QAM_Error_Symbol)],QAM_Error_Symbol,'b*-.','MarkerSize',10);
hold on;
semilogy([1:length(QAM_Error_Symbol_theory)],QAM_Error_Symbol_theory,'rO:','MarkerSize',10);
title('16-QAM SNR to Symbol Error Rate');
legend('Simulated','Theoretical');
grid on
hold off;