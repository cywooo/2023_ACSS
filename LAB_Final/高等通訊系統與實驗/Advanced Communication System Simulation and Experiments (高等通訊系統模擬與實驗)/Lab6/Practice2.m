clear all;
close all;
%% Generate a sinusoidal signal, downsample the signal, and then upsample the downsample signal
% Generate a sinusoidal signal
N = 128;
f = 1/10;
s = cos(2*pi*f*[1:N]);
figure(1);
stem(s);
title('Sinusoidal signal');
figure(2);
plot(fftshift(abs(fft(s))));
title('Spectrum of Sinusoidal signal');

% Downsample the signal
M = 3;
downsample_s = s(1:M:end);
figure(3);
stem(downsample_s);
title('Downsample sinusoidal signal');
figure(4);
plot(fftshift(abs(fft(downsample_s))));
title('Downsample sinusoidal spectrum');

% Upsample the downsample signal
L = 3;
upsample_s = zeros(1,N);
upsample_s(1:L:end) = downsample_s;
figure(5);
stem(upsample_s);
title('Upsample sinusoidal signal');
figure(6);
plot(fftshift(abs(fft(upsample_s))));
title('Upsample sinusoidal spectrum');

%% Design a an FIR LPF and let the upsampled signal pass the filter such that the upsampled signal is similar to the original signal
% FIR
recovered_s_FIR = filter(FIR_LPF,upsample_s);
figure(7);
stem(L*recovered_s_FIR(51:end)); % Look at the delay
figure(8);
plot(fftshift(abs(fft(recovered_s_FIR))));

% IIR
recovered_s_IIR = filter(IIR_LPF,upsample_s);
figure(9);
stem(L*recovered_s_IIR(40:end));
figure(10);
plot(fftshift(abs(fft(recovered_s_IIR))));
