clear all
close all

%% Transmitter
% BPSK
N = 1024;
symbol = randi([0 1],1,N);
symbol(symbol==0) = -1;

% Upsample
Lu = 16;
upsample_sym = zeros(1,Lu*length(symbol));
upsample_sym(1:Lu:end) = symbol;

% Gaussian filter
M = Lu;
BT = 0.5;
trun = 5;
gaussian_sym = conv(upsample_sym,Gaussian_filter(M,BT,trun),'same');

% Sum
fd = 150*10^3;
fb = 1*10^6;
fs = fb*Lu;
Tb = 1/(fb*Lu);
sum_sf = zeros(1,length(gaussian_sym));
sum_sf(1) = gaussian_sym(1);
for m=2:length(gaussian_sym)
    sum_sf(m) = sum_sf(m-1) + gaussian_sym(m);
end
phase_sf = exp(i*2*pi*fd*Tb.*sum_sf);

% IF_upconvert
f_IF = 2*10^6;
t = 0:length(phase_sf)-1;
w_IF = 2*pi*f_IF/fs;
IF_sym = real(phase_sf.*exp(i*w_IF*t));

% DAC
Ldac = 4;
dac_sym = zeros(1,Ldac*length(IF_sym));
dac_sym(1:Ldac:end) = IF_sym;

% DMA
dma_sym = filter(DMA_new,dac_sym);
gd = 10;
dma_sym = dma_sym(gd+1:end);

% % Up-convert
% fc = 8*10^6;
% wc = 
% x = dma_sym *exp(
%% Plot
% xaxis = (1/64/10^6) .* [f_IF-3*10^6, f_IF-2.5*10^6, f_IF-1.5*10^6, f_IF-1*10^6, f_IF, f_IF+1*10^6, f_IF+1.5*10^6, f_IF+2.5*10^6, f_IF+3*10^6, f_IF+4*10^6];
% yaxis = [-66,-46, -26, 0, 0, 0, -26, -46, -66,-66];

[px,f] = pwelch(dma_sym,[],[],[],1);

p1 = (f_IF-2.5*10^6)/64/10^6;
p2 = (f_IF-1.5*10^6)/64/10^6;
p3 = (f_IF-1*10^6)/64/10^6;
p4 = (f_IF+1*10^6)/64/10^6;
p5 = (f_IF+1.5*10^6)/64/10^6;
p6 = (f_IF+2.5*10^6)/64/10^6;

f0 = find(f<p1);
f1 = find(f>=p1 & f<p2);
f2 = find(f>=p2 & f<p3);
f3 = find(f>=p3 & f<p4);
f4 = find(f>=p4 & f<p5);
f5 = find(f>=p5 & f<p6);
f6 = find(f>=p6);

mask(f0) = max(10*log10(px))-66;
mask(f1) = max(10*log10(px))-46;
mask(f2) = max(10*log10(px))-26;
mask(f3) = max(10*log10(px));
mask(f4) = max(10*log10(px))-26;
mask(f5) = max(10*log10(px))-46;
mask(f6) = max(10*log10(px))-66;

plot(f,10*log10(px));
hold on
plot(f,mask);