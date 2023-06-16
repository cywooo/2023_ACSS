clear all;
close all;
%% Create an arbitrary digital signal
N = -64:63;
w = 64;
tri = tripuls(N,w);
figure(1);
stem(tri);
title('Triangular pulse');
figure(2);
plot(fftshift(abs(fft(tri))));
title('Spectrum of Trianglular pulse');

%% Design a DMA filter and let it pass through an ideal DAC with the upsampling rate of 32.
% Upsampling(DAC)
L = 32;
upsample_length = L * length(tri);
upsample_tri = zeros(1,upsample_length);
upsample_tri(1:L:end) = tri;
figure(3);
stem(upsample_tri);
title('Upsample Trianglular pulse');
figure(4);
plot(fftshift(abs(fft(upsample_tri))));
title('Spectrum of Upsample Trianglular pulse');

% DMA filtering
analog_tri = filter(DMA,upsample_tri);
figure(5);
plot(L*analog_tri);
title('Analog Trianglar pulse');
figure(6);
plot(fftshift(abs(fft(analog_tri))));
title('Spectrum of Analog Trianglar pulse');

%% Let the modeled analog signal pass with an ADC with the rate of 32
M = 32;
downsample_tri = L*analog_tri(1:M:end);
figure(7);
stem(downsample_tri);
title('Downsampled Analog signal');
figure(8);
plot(fftshift(abs(fft(downsample_tri))));
title('Spectrum of Downsampled Analog signal');