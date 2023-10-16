clear;close all;clc;
load('P1');

% -----------------------------(a)---------------------------------
Es = mean((abs(g1a).^2));
En = Es/10;
v = sqrt(En)*randn(1,177);
r1a_np = g1a+v;

% -----------------------------(b)---------------------------------
% FIR filter
rlb_np_f = fft(r1a_np);
r1b_fp = filter(r1b_filter,r1a_np);
r1b_fp_hat = filter(r1b_filter,g1a);
% snr
r1b_snr = 10*log10(mean(abs(r1b_fp).^2) / mean(abs(r1b_fp-r1b_fp_hat).^2));

% % plot sig before filtering
% f = linspace(-1/2,1/2,length(rlb_np_f));
% figure();plot(f,fftshift(abs(rlb_np_f)));
% plot sig after filtering
figure();plot(g1a(1:end-11));hold on;plot(r1b_fp(11:end));

% -----------------------------(c)---------------------------------
rlc_np_f = fft(r1a_np);
r1c_fp = filter(r1c_filter,r1a_np);
r1c_fp_hat = filter(r1c_filter,g1a);
% snr
r1c_snr = 10*log10(mean(abs(r1c_fp).^2) / mean(abs(r1c_fp-r1c_fp_hat).^2));

% % plot sig before filtering
% f = linspace(-1/2,1/2,length(rlc_np_f));
% figure();plot(f,fftshift(abs(rlc_np_f)));
% plot sig after filtering
figure();plot(g1a(1:end-30));hold on;plot(r1c_fp(30:end));