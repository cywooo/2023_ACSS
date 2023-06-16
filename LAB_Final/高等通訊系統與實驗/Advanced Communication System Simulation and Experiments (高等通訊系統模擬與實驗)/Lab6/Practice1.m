clear all;
close all;
%% Generate a sinusoidal signal, downsample the signal, and observe its spectrum
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

% Downsample and observe its spectrum
M = 3;
downsample_s = s(1:M:end);
figure(3);
stem(downsample_s);
title('Downsample sinusoidal signal');
figure(4);
plot(fftshift(abs(fft(downsample_s))));
title('Downsample sinusoidal spectrum');

%% Upsample the downsample signal, and observe its spectrum
L = 3;
upsample_s = zeros(1,N);
upsample_s(1:L:end) = downsample_s;
figure(5);
stem(upsample_s);
title('Upsample sinusoidal signal');
figure(6);
plot(fftshift(abs(fft(upsample_s))));
title('Upsample sinusoidal spectrum');

