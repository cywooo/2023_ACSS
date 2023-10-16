clear all
close all

%% Input signal
N = 1024;
data = randi([0 1],1,N);
data = reshape(data,2,N/2);
symbol = exp((1i*pi/4) .* ((-1).^(data(1,:))) .* (3.^(data(2,:))));

%% Upsampling
Lu = 16;
upsample_sym = zeros(1,length(symbol)*Lu);
upsample_sym(1:Lu:end) = symbol;

%% Gaussian filter
M = Lu;
BT = 0.5;
trun = 5;
gaussian_sym = conv(upsample_sym,Gaussian_filter(M,BT,trun),'same');

%% Summation to phase (Integration)
Sf = zeros(1,length(gaussian_sym));
for m=1:length(gaussian_sym)
    Sf(m) = sum(gaussian_sym(1:m));
end
fd = 150*1e3;
fb = 1*1e6;
Tb = 1/(fb*Lu);
Sf_r = real(Sf);
Sf_i = imag(Sf);
%integrate_sym = exp(1i*2*pi*fd*Tb.*Sf);
integrate_sym_r = exp(1i*2*pi*fd*Tb.*Sf_r);
integrate_sym_i = exp(1i*2*pi*fd*Tb.*Sf_i);

%% IF modulation
f_IF = 2*1e6;
fs = fb*Lu;
m = 0:length(integrate_sym_r)-1;
IF_sym_r = real(integrate_sym_r .* exp(1i*2*pi*f_IF/fs*m));
IF_sym_i = real(integrate_sym_i .* exp(1i*2*pi*f_IF/fs*m));

%% DAC
Ldac = 4;
dac_sym_r = zeros(1,Ldac*length(IF_sym_r));
dac_sym_i = zeros(1,Ldac*length(IF_sym_i));
dac_sym_r(1:Ldac:end) = IF_sym_r;
dac_sym_i(1:Ldac:end) = IF_sym_i;

%% DMA
dma_sym_r = filter(DMA_new,dac_sym_r);
group_delay = 22;
dma_sym_r = dma_sym_r(group_delay+1:end);
dma_sym_i = filter(DMA_new,dac_sym_i);
dma_sym_i = dma_sym_i(group_delay+1:end);

%% DMA (Receiver)
dma_received_r = filter(DMA_new,dma_sym_r);
dma_sym_r = dma_sym_r(group_delay+1:end);
dma_received_i = filter(DMA_new,dma_sym_i);
dma_sym_i = dma_sym_i(group_delay+1:end);

%% ADC
Ladc = Ldac;
adc_sym_r = dma_received_r(1:Ladc:end);
adc_sym_i = dma_received_i(1:Ladc:end);

%% IF demodulation
m = 0:length(adc_sym_r)-1;
IF_received_r = adc_sym_r .* exp(-1i*2*pi*f_IF/fs*m);
m = 0:length(adc_sym_i)-1;
IF_received_i = adc_sym_i .* exp(-1i*2*pi*f_IF/fs*m);

%% SRRC filter
trun = 5;
M = Lu;
alpha = 0.3;
srrc_sym_r = conv(IF_received_r,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
srrc_sym_r = srrc_sym_r(delay+1:end-delay);
srrc_sym_i = conv(IF_received_i,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
srrc_sym_i = srrc_sym_i(delay+1:end-delay);

%% Phase
phase_sym_r = unwrap(angle(srrc_sym_r));%/(2*pi*fd*Tb);
phase_sym_i = unwrap(angle(srrc_sym_i));%/(2*pi*fd*Tb);
phase_sym = (phase_sym_r + 1i*phase_sym_i)/(2*pi*fd*Tb);
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
downsample_received = gaussian_received(1:Ld:end)*1.25;

%% Detection
detection = (1/sqrt(2)) * (((real(downsample_received)>0)-(real(downsample_received)<0)) + 1i*((imag(downsample_received)>0)-(imag(downsample_received)<0)));

%% EVM
EVM = 10*log10(mean(abs(symbol-downsample_received).^2) / mean(abs(symbol).^2));

%% Plot
figure(1)
stem(real(symbol));
hold on
stem(real(downsample_received));

figure(2);
plot(Sf_r);
hold on
plot(real(phase_sym_r)/(2*pi*fd*Tb));

% figure(3);
% [px,f] = pwelch(dma_sym,[],[],[],1);
% plot(f,10*log10(px));
% 
% figure(4);
% p1 = (f_IF-2.5*1e6)/(1*1e6*Lu*Ldac);
% p2 = (f_IF-1.5*1e6)/(1*1e6*Lu*Ldac);
% p3 = (f_IF-1*1e6)/(1*1e6*Lu*Ldac);
% %p4 = (f_IF)/(1*1e6*Lu*Ldac);
% p4 = (f_IF+1*1e6)/(1*1e6*Lu*Ldac);
% p5 = (f_IF+1.5*1e6)/(1*1e6*Lu*Ldac);
% p6 = (f_IF+2.5*1e6)/(1*1e6*Lu*Ldac);
% 
% f0 = find(f<p1);
% f1 = find(f>=p1 & f<p2);
% f2 = find(f>=p2 & f<p3);
% f3 = find(f>=p3 & f<p4);
% f4 = find(f>=p4 & f<p5);
% f5 = find(f>=p5 & f<p6);
% f6 = find(f>=p6);
% 
% mask(f0) = max(10*log10(px))-66;
% mask(f1) = max(10*log10(px))-46;
% mask(f2) = max(10*log10(px))-26;
% mask(f3) = max(10*log10(px));
% mask(f4) = max(10*log10(px))-26;
% mask(f5) = max(10*log10(px))-46;
% mask(f6) = max(10*log10(px))-66;
% 
% plot(f,10*log10(px));
% hold on
% plot(f,mask);
% axis([0 0.12 -80 40]);



