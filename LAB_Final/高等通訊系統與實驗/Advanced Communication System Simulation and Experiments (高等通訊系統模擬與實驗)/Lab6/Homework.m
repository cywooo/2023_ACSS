clear all;
close all;
%% For a given signal, try to conduct downsampling with a largest factor without causing distortion
% Load given signal
A1 = load('A1.mat');
a1 = A1.a1;
N = length(a1);
fs = 1;
f = 2*pi/pi*fs*((-N/2):(N/2)-1)/N; % pi
figure(1);
plot(a1);
title('Input Signal');
figure(2);
plot(f,fftshift(abs(fft(a1))));
title('Spectrum of Input Signal');

% Downsampling with a largest factor
M = 3;
downsample_signal = a1(1:M:end);
figure(3);
plot(downsample_signal);
title('Downsampled Signal');
Nd = length(downsample_signal);
fd = fs*((-Nd/2):(Nd/2)-1)/Nd; % pi
figure(4);
plot(fd,fftshift(abs(fft(downsample_signal))));
title('Spectrum of Downsampled Signal');

%% Design a filter and conduct upsampling to recover the signal
% Upsampling
L = M;
upsample_length = L * length(downsample_signal);
upsample_signal = zeros(upsample_length,1);
upsample_signal(1:L:end) = downsample_signal;
figure(5);
plot(upsample_signal);
title('Upsampled Signal');
Nu = length(upsample_signal);
fu = 2*pi/pi*fs*((-Nu/2):(Nu/2)-1)/Nu; % pi
figure(6);
plot(fu,fftshift(abs(fft(upsample_signal))));
title('Spectrum of Upsampled Signal');

% Filtering
recovered_signal = filter(IIR_LPF_hw,upsample_signal);
figure(7);
plot(L*recovered_signal(10:end));
title('Recovered Signal');
Nr = length(recovered_signal);
fr = 2*pi/pi*fs*((-Nr/2):(Nr/2)-1)/Nr; % pi
figure(8);
plot(fr,fftshift(abs(fft(L*recovered_signal))));
title('Spectrum of Recovered Signal');
figure(9);
plot(a1);
hold on
plot(L*recovered_signal(10:end));
legend('origin','recover');
title('Compare of Origin and Recovered Signal');