clear all; clc; close all;

SNR = [-10:5:20]; % inf;%    
batch = 500;
%% 
for q = [1:length(SNR)]
    for qq = [1:batch]
        signal_length = 1e4; % signal length adjust
        % signal generate------
        s_2 = sign(randn(1,signal_length)); %
        
        x = s_2;
        % Channel===============                        
        h_2_to_2 = [ randn(1,1)+randn(1,1)*1i randn(1,1)+randn(1,1)*1i ;...
                     randn(1,1)+randn(1,1)*1i randn(1,1)+randn(1,1)*1i ];           

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

        s_11_ZF_reDig = mean(y_ZF);

        s_11_MMSE_reDig = mean(y_MMSE);

        OSNR_ZF = 1/(sum(abs(real(s_11_ZF_reDig)-s_2).^2)/(signal_length));
        OSNR_ZF_dB(q,qq) = 10*log10(OSNR_ZF);
        OSNR_MMSE = 1/(sum(abs(real(s_11_MMSE_reDig)-s_2).^2)/(signal_length));
        OSNR_MMSE_dB(q,qq) = 10*log10(OSNR_MMSE);

        BER_ZF(q,qq) = sum(abs(sign(real(s_11_ZF_reDig))-s_2))/2/(signal_length);
        BER_MMSE(q,qq) = sum(abs(sign(real(s_11_MMSE_reDig))-s_2))/2/(signal_length);
    end
    
end
OSNR_ZF_dB_ = sum(OSNR_ZF_dB')'/batch;
OSNR_MMSE_dB_ = sum(OSNR_MMSE_dB')'/batch;
BER_ZF_ = sum(BER_ZF')'/batch;
BER_MMSE_ = sum(BER_MMSE')'/batch;

figure(1);
semilogy(SNR,BER_ZF_,"b");
hold on;
semilogy(SNR,BER_MMSE_,"r");
hold off;
legend("ZF","MMSE");
xlabel("SNR");
ylabel("BER");

figure(2);
plot(SNR,OSNR_ZF_dB_,"b");
hold on;
plot(SNR,OSNR_MMSE_dB_,"r");
hold off;
legend("ZF","MMSE");
xlabel("SNR");
ylabel("SNR_o");