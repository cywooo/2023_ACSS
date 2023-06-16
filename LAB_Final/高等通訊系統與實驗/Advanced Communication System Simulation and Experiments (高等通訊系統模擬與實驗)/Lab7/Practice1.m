clear all;
close all;

alpha = 0.3;
M = 8;
trun = 5; % See how many T we want to generate
n = [-M*trun:M*trun]+0.00001;
% RC
RC = sinc(n/M) .* cos(2*pi*alpha*n/2/M) ./ (1 - 16*alpha^2*n.^2/4/M^2);

% SRRC
nom = (4*alpha/pi) * (cos((1+alpha)*pi*n/M) + M*sin((1-alpha)*pi*n/M)/4/alpha./n);
%nom = (4*alpha/pi) * (cos((1+alpha)*pi*n/M) + sinc((1-alpha)*n/M)/4/alpha*(1-alpha)*pi);
den = 1-(4*alpha*n/M).^2;
SRRC = nom./den;

% Convolution
conv_RC = conv(SRRC,SRRC,'same'); % conv_RC = conv(RC,RC);
y = conv_RC/max(conv_RC);
%practical_SRRC = conv(ones(1,M),conv_RC);

figure(1);
subplot(3,1,1);
stem(RC);
title('RC pulse');
subplot(3,1,2);
stem(SRRC_filter(trun,M,alpha));
title('SRRC pulse');
subplot(3,1,3);
stem(y);
title('Convolution of SRRC pulse');

figure(2);
subplot(3,1,1);
stem(fftshift(abs(fft(RC))));
title('Spectrum of RC pulse');
subplot(3,1,2);
stem(fftshift(abs(fft(SRRC))));
title('Spectrum of SRRC pulse');
subplot(3,1,3);
stem(fftshift(abs(fft(y))));
title('Spectrum of Convolution of SRRC pulse');

figure(3);
stem(fftshift(abs(fft(RC))));
hold on
stem(fftshift(abs(fft(SRRC))));
stem(fftshift(abs(fft(y))));
legend('RC','SRRC','Conv');

