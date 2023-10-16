clear all
close all

%% BPSK
N = 64;
symbol = randi([0 1],2,N);
symbol(symbol==0) = -1;
symbol_1 = symbol(1,:);
symbol_2 = symbol(2,:);

%% Channel
nt = 2;
nr = 2;
H = 1/sqrt(2)*(randn(nr,nt) + 1i*randn(nr,nt));
H = H ./ sqrt(sum(abs(H).^2));
channel_sym = H*symbol;

%% Noise
SNR = 10;
signal_power = mean(abs(symbol_1).^2);
noise_power = signal_power*10^(-SNR/10);
noise = sqrt(noise_power/2) * (randn(2,length(channel_sym)) + 1i*randn(2,length(channel_sym)));
received_sym = channel_sym + noise;

%% ZF detection
ZF_detector = inv(H);
ZF_recover = real(ZF_detector*received_sym);

%% MMSE detection
rho = [signal_power/noise_power signal_power/noise_power];
MMSE_detector = inv(H'*H + inv(diag(rho)))*H';
MMSE_recover = real(MMSE_detector*received_sym);

%% SNR
symbol_power = mean(abs(symbol).^2,'all');
error_power_ZF = mean(abs(ZF_recover-symbol).^2,'all');
ZF_SNR = 10*log10(symbol_power/error_power_ZF)

error_power_MMSE = mean(abs(MMSE_recover-symbol).^2,'all');
MMSE_SNR = 10*log10(symbol_power/error_power_MMSE)

%% Plot
figure(1);
stem(symbol_1);
hold on
stem(ZF_recover(1,:));
title('ZF recover of first signal');
figure(2);
stem(symbol_2);
hold on
stem(ZF_recover(2,:));
title('ZF recover of second signal');

figure(3);
stem(symbol_1);
hold on
stem(MMSE_recover(1,:));
title('MMSE recover of first signal');
figure(4);
stem(symbol_2);
hold on
stem(MMSE_recover(2,:));
title('MMSE recover of second signal');
