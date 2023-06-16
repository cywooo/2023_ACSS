clear all
close all

%% Input signal
N = 1024;
symbol = randi([0 1],1,N);
symbol(symbol==0) = -1;

%% Upsampling
Lu = 16;
upsample_sym = zeros(1,length(symbol)*Lu);
upsample_sym(1:Lu:end) = symbol;

%% Gaussian filter
M = Lu;
BT = 0.5;
trun = 5;
gaussian_sym = conv(upsample_sym,Gaussian_filter(M,BT,trun),'same');
% delay = (length(Gaussian_filter(M,BT,trun))-1)/2;
% gaussian_sym = gaussian_sym(delay+1:end-delay);

%% Summation to phase (Integration)
Sf = zeros(1,length(gaussian_sym));
for m=1:length(gaussian_sym)
    Sf(m) = sum(gaussian_sym(1:m));
end
fd = 150*1e3;
fb = 1*1e6;
Tb = 1/(fb*Lu);
integrate_sym = exp(1i*2*pi*fd*Tb.*Sf);

%% IF modulation
f_IF = 2*1e6;
fs = fb*Lu;
m = 0:length(integrate_sym)-1;
IF_sym = real(integrate_sym .* exp(1i*2*pi*f_IF/fs*m));

%% DAC
Ldac = 4;
dac_sym = zeros(1,Ldac*length(IF_sym));
dac_sym(1:Ldac:end) = IF_sym;

%% DMA
trun = 5;
M = Ldac;
alpha = 0.1;
dma_sym = conv(dac_sym,SRRC_filter(trun,M,alpha),'same');
% delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
% dma_sym = dma_sym(delay+1:end-delay);

%% DMA (Receiver)
dma_received = conv(dma_sym,SRRC_filter(trun,M,alpha),'same');
% delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
% dma_received = dma_received(delay+1:end-delay);

%% ADC
Ladc = Ldac;
adc_sym = dma_received(1:Ladc:end);

%% IF demodulation
m = 0:length(adc_sym)-1;
IF_received = adc_sym .* exp(-1i*2*pi*f_IF/fs*m);

%% SRRC filter
M = Lu;
srrc_sym = conv(IF_received,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
srrc_sym = srrc_sym(delay+1:end-delay);

%% Phase
phase_sym = unwrap(angle(srrc_sym))/(2*pi*fd*Tb);
phase_received = zeros(1,length(phase_sym));
phase_received(1) = phase_sym(1);
for m=2:length(phase_sym)
    phase_received(m) = phase_sym(m) - phase_sym(m-1);
end

%% Gaussian filter
M = Lu;
gaussian_received = conv(phase_received,Gaussian_filter(M,BT,trun),'same');

%% Downsampling
Ld = Lu;
downsample_received = gaussian_received(1:Ld:end);

%% Detection
detection = (downsample_received>0) - (downsample_received<0);

%% Plot
figure(1)
stem(symbol);
hold on
stem(detection);

figure(2);
plot(Sf);
hold on
plot(phase_sym);

figure(3);
[px,f] = pwelch(dma_sym,[],[],[],1);
plot(f,10*log10(px));

figure(4);
p1 = (f_IF-2.5*1e6)/(1*1e6*Lu*Ldac);
p2 = (f_IF-1.5*1e6)/(1*1e6*Lu*Ldac);
p3 = (f_IF-1*1e6)/(1*1e6*Lu*Ldac);
%p4 = (f_IF)/(1*1e6*Lu*Ldac);
p4 = (f_IF+1*1e6)/(1*1e6*Lu*Ldac);
p5 = (f_IF+1.5*1e6)/(1*1e6*Lu*Ldac);
p6 = (f_IF+2.5*1e6)/(1*1e6*Lu*Ldac);

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
axis([0 0.12 -80 40]);



