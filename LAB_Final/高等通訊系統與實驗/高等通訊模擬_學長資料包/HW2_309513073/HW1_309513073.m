
f = 1000000;
fs = 4000000;
N = 256;  %DFT size N
t = [0:N-1];

x=cos(2*pi*f/fs*t);   %0.25M=250000

% figure(1)
% plot(x)
xf =abs(fft(x));
figure(1)
plot(abs(fft(x)))

fs2 = 1500000;      %1/1.5
x2=cos(2*pi*f/fs2*t);

% figure(1)
% stem(x2)
figure(2)
plot(abs(fft(x2)))





