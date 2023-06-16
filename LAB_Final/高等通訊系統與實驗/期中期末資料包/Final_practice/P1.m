clear;close all;
rng(223);

Nt = 2;
QPSK_seq = [-1-1j -1+1j 1-1j 1+1j];
QPSK = QPSK_seq(randi(4,[Nt,64]));
seq_idx = randi(4,[Nt,64]); 
seq_idx_1 = seq_idx(1,:);
seq_idx_2 = seq_idx(2,:);

% QPSK = randsample(QPSK_seq,64,true);
QPSK_1 = QPSK(1,:);
QPSK_2 = QPSK(2,:);

%----------------------------------------channel
Nt = 2;
Nr = 2;
H = 1/sqrt(2)*(randn(Nr,Nt)+1i*randn(Nr,Nt));
H = H./sqrt(sum(abs(H).^2));
receive_signal = H * QPSK;

%----------------------------------------add noise
SNR_dB = 10;
signal_power = mean(abs(QPSK_1).^2);
N0 = signal_power*10^(-SNR_dB/10);
noise = sqrt(N0/2)*(randn(Nr,length(receive_signal))+1i*randn(Nr,length(receive_signal)));
receive_signal = noise + receive_signal;

%-----------------------------------------signal detection
ZF = inv(H);
rho = [signal_power/N0 signal_power/N0];
MMSE = inv(H'*H+diag(rho))*H';

detect_signal_ZF = ZF * receive_signal;
detect_signal_MMSE = MMSE * receive_signal;

[noneed1,idx_ZF1] = min(abs(detect_signal_ZF(1,:) - QPSK_1.'));
symbol_error_ZF1 = sum(idx_ZF1 ~= seq_idx_1)/64;
[noneed2,idx_ZF2] = min(abs(detect_signal_ZF(2,:) - QPSK_2.'));
symbol_error_ZF2 = sum(idx_ZF2 ~= seq_idx_2)/64;
symbol_error_ZF = symbol_error_ZF1 + symbol_error_ZF2;



