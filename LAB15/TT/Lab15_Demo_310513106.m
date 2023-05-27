clear all;close all;
%% setting
fs = 1e6;fd = 150e3;
f_IF = 2e6;f_adc = 16e6;f_dac = f_adc;
fan = 64e6;

%gaussian filter
BT =.5;
M = f_dac/fs;
T = 1;Tb = M*T;span = 7;
t_gau = linspace(-16, span*T, span*2*1*T+1);
gau = exp(-2*pi^2 / log10(2) * (1/2/M)^2 .* t_gau.^2);
gau = gau * (BT/Tb)*sqrt(2*pi/log10(2));

figure;
stem(t_gau,gau);

N = 100;
a = 2*round(rand(1, N)) - 1;
aa = zeros(1, N*M);
aa(1:M:end) = a;
take = conv(gau, aa);
a_gau = take(floor(length(gau)/2)+1 : end - floor(length(gau)/2));
%% sum
a_sum = zeros(1, length(a_gau));
tmpp = 0;
for i = 1:length(a_gau)
    tmpp = a_gau(i) + tmpp;
    a_sum(i) = tmpp;
end

modu_in = fd *(1/fs)/M;
a_exp = exp(1j*2*pi*modu_in*a_sum);
%% IF UP
t_IF = 1:length(a_exp);
carr_IF = exp(1j*2*pi*f_IF/f_dac*t_IF);
a_IF = real(a_exp .* carr_IF);
figure;
plot([1:length(a_IF)],a_IF);
%% DAC
u_dac = fan / f_dac;
au_dac = zeros(1, length(a_IF)*u_dac);
au_dac(1:u_dac:end) = a_IF;

figure, plot(linspace(-1, 1, length(au_dac)), abs(fftshift(fft(au_dac))));

a_dac = filter(DMA, au_dac); 
%%
[px, f] = pwelch(a_dac, [], [], [], 1);
figure;
plot(f, 10*log10(px));
hold on ;
mask = (-40).*(f<=(2/64-2.5/64)) + (-20).*(f<=(2/64-1.5/64)).*(f>(2/64-2.5/64)) + 0.*(f<=(2/64-1/64)).*(f>(2/64-1.5/64)) + 26.*(f<=(2/64+1/64)).*(f>(2/64-1/64)) + 0.*(f<=(2/64+1.5/64)).*(f>(2/64+1/64)) + (-20).*(f<=(2/64+2.5/64)).*(f>(2/64+1.5/64)) + (-40).*(f>(2/64+2.5/64));
plot(f, mask);
hold off;