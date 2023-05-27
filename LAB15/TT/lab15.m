clc 
clear all

%% setting
fs = 1e6;
fd = 150e3;
f_IF = 2e6;
f_adc = 16e6;
f_dac = f_adc;
fan = 64e6;

%% Gaussian filter
BT = .5;
M = f_dac / fs;
Tg = 1;
B = BT / Tg;
C = sqrt(2*pi/log10(2)) * B;
Tb = M*Tg;
span = 7;
t_gau = linspace(-span*Tg, span*Tg, span*2*1*Tg+1);
gau = C * exp(-2*pi^2 / log(2) * (1/2/M)^2 .* t_gau.^2);
gau = gau / norm(gau);
%
% figure, plot(gau), title('Gaussian filter')

%% Origial
N = 6e2;
a = 2*round(rand(1, N)) - 1;
aa = zeros(1, N*M);
aa(1:M:end) = a;
%
% figure(1), plot(linspace(-1, 1, length(aa)), abs(fftshift(fft(aa))))

take = conv(gau, aa);
a_gau = take(floor(length(gau)/2)+1 : end - floor(length(gau)/2));
%
% figure(2), plot(linspace(-1, 1, length(a_gau)), abs(fftshift(fft(a_gau))))

%% Sum
a_sum = zeros(1, length(a_gau));
tmpp = 0;
for i = 1:length(a_gau)
    tmpp = a_gau(i) + tmpp;
    a_sum(i) = tmpp;
end

%% Phase it
% modu_in = fd / f_dac * Tb; % modulation index
modu_in = fd *(1/fs)/M;
% modu_in = .6;
a_exp = exp(1j*2*pi*modu_in*a_sum);

%% IF_up
t_IF = 1:length(a_exp);
carr_IF = exp(1j*2*pi*f_IF/f_dac*t_IF);
a_IF = real(a_exp .* carr_IF);
% a_IF = (a_exp .* carr_IF);
%
% figure, plot(linspace(-1, 1, length(a_IF)), abs(fftshift(fft(a_IF)))),
% title('Before DAC filter')

%% DAC
u_dac = fan / f_dac;
au_dac = zeros(1, length(a_IF)*u_dac);
au_dac(1:u_dac:end) = a_IF;

% figure, plot(linspace(-1, 1, length(au_dac)), abs(fftshift(fft(au_dac))))
% title('After DAC filter')

a_dac = filter(dma, au_dac); % d = 4

% figure, plot(linspace(-1, 1, length(a_dac)), abs(fftshift(fft(a_dac))))
% title('Before up-conversion in frequency domain')

[px, f] = pwelch(a_dac, [], [], [], 1);
figure
plot(f, 10*log10(px))
hold on 
mask = (-40).*(f<=(2/64-2.5/64)) + (-20).*(f<=(2/64-1.5/64)).*(f>(2/64-2.5/64)) + 0.*(f<=(2/64-1/64)).*(f>(2/64-1.5/64)) + 26.*(f<=(2/64+1/64)).*(f>(2/64-1/64)) + 0.*(f<=(2/64+1.5/64)).*(f>(2/64+1/64)) + (-20).*(f<=(2/64+2.5/64)).*(f>(2/64+1.5/64)) + (-40).*(f>(2/64+2.5/64));

plot(f, mask)
hold off
% mask = (0:1/128:1/2);
% mask
%{

%% Up-conversion
fc = 32e6;
fn = (fc - f_IF) / (M*u_dac*fs);
tu = 1:length(a_dac);
carr_up = exp(1j*2*pi*fn*tu);
a_up = a_dac .* carr_up;
%
% figure, plot(linspace(-1, 1, length(a_up)), abs(fftshift(fft(a_up))))
% title('After up-conversion in frequency domain')

[px, f] = pwelch(a_up, [], [], [], 1);
figure
plot(f, 10*log10(px))


%% Noise
SNR = 15;
snr = 10^(SNR/10);
noise = randn(1, length(a_up)) * sqrt(mean(a.^2) / snr);
a_noise = a_up;% + noise;
%
% figure, plot(linspace(-1, 1, length(a_noise)), abs(fftshift(fft(a_noise))))
% title('Add noise with 15dB in frequency domain')


%% Down-conversion
carr_down = exp(-1j*2*pi*fn*tu);
a_down = a_noise .* carr_down;


%% ADC
a_adc = filter(dma, a_down); % d = 4

%% down-sampling
ad = a_adc(8:u_dac:end);

%% IF_down
t_IF_d = 1:length(ad);
carr_IF_d = exp(-1j*2*pi*f_IF/f_adc*t_IF_d);
ad_IF = ad .* carr_IF_d;

%% sccr
n = u_dac; % 4
alpha = .1;
T = 1;
span_srrc = 8;
t = linspace(-span_srrc*T, span_srrc*T, 2*span_srrc*T*n+1) - 0.000001;
srrc = 4 * alpha/ pi * (cos((1+alpha)* pi* t/ T) + T* sin((1-alpha)* pi* t/ T)./ (4* alpha* t))./ (1- (4*alpha*t/T).^ 2);
srrc = srrc/norm(srrc);
%
% figure(5), plot(linspace(-1, 1, length(srrc)), abs(fftshift(fft(srrc))))
%% SRRC filter
a_srrc = conv(ad_IF, srrc);
a_srrc = a_srrc(floor(length(srrc)/2)+1 : end - floor(length(srrc)/2));
a_phase = unwrap(angle(a_srrc));
%
% figure(3), plot(linspace(-1, 1, length(ad_IF)), abs(fftshift(fft(ad_IF))))
% figure(4), plot(linspace(-1, 1, length(a_srrc)), abs(fftshift(fft(a_srrc))))

%% Diff
% a_diff = [a_srrc(1) diff(a_srrc)];
a_diff = diff(a_phase);
a_diff = [a_diff 0];

%% Gaussian filter
ad_gau = conv(gau, a_diff);
ad_gau = ad_gau(floor(length(gau)/2)+1 : end - floor(length(gau)/2));

%% Down-sampling
s = ad_gau(1:M:end);

% figure
% stem(a)
% hold on
% % stem(s)
% stem(sign(s), 'x')
% hold off
% legend('original signal', 'recovered signal')

figure, plot(a_sum), hold on, plot(a_phase), hold off, legend('modulated', 'demodulated'), title('Decrease modulation index(h = 0.1)')

error_rate = sum(sign(s) ~= a)/N

%}