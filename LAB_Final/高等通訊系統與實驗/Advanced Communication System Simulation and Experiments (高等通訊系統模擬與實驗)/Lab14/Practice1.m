clear all
close all

%% SQNR of white uniform signal
N = 1024;
x = 2*rand(1,N)-1;
numOfBits = 6;
dynamic_range = 1;

xq = zeros(1,N);
for n=1:N
    xq(n) = ADC(x(n),numOfBits,dynamic_range);
end

SQNR_1 = 10*log10(mean(abs(x).^2)/mean(abs(x-xq).^2));
disp(['SQNR of white uniform signal is ' ,num2str(SQNR_1)]);

%% 6dB Rule
numOfBits = 5;

xq = zeros(1,N);
for n=1:N
    xq(n) = ADC(x(n),numOfBits,dynamic_range);
end

SQNR_2 = 10*log10(mean(abs(x).^2)/mean(abs(x-xq).^2));
SQNR = SQNR_1-SQNR_2;
disp(['Verify of 6dB rule: ' ,num2str(SQNR)]);

%% SQNR of white Gaussian signal
x = randn(1,N);% + 1i*randn(1,N);
numOfBits = 4;
dynamic_range = 8;
y = linspace(-dynamic_range,dynamic_range,N);

xq = zeros(1,N);
for n=1:N
    xq(n) = ADC(x(n),numOfBits,dynamic_range);
    yq(n) = ADC(y(n),numOfBits,dynamic_range);
end

SQNR = 10*log10(mean(abs(x).^2)/mean(abs(x-xq).^2));
disp(['SQNR of white Gaussian signal is ' ,num2str(SQNR)]);

figure
plot(y,yq);
