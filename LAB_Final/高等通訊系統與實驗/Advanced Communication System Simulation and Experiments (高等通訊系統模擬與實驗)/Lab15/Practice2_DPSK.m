clear all
close all

%% Input signal
N = 1024;
input_bit = randi([0 1],1,N);
data = reshape(input_bit,2,N/2);

%% Mapping
symbol = (pi/4) .* ((-1).^(data(1,:))) .* (3.^(data(2,:)));

%% Phase Integration
Sf = zeros(1,length(symbol));
for m = 1:length(symbol)
    Sf(m) = sum(symbol(1:m));
end
integrate_phase = exp(1i*Sf);

%% Upsampling
Lu = 16;
upsample_sym = zeros(1,length(integrate_phase)*Lu);
upsample_sym(1:Lu:end) = integrate_phase;

%% Pulse shaping
trun = 5;
M = Lu;
alpha = 0.4;
srrc_sym = conv(upsample_sym,SRRC_filter(trun,M,alpha),'same');

%% IF modulation
f_IF = 2*1e6;
fb = 1*1e6;
fs = fb*Lu;
m = 0:length(srrc_sym)-1;
IF_sym = srrc_sym .* exp(1i*2*pi*f_IF/fs*m);

%% DAC
Ldac = 4;
dac_sym = zeros(1,length(IF_sym)*Ldac);
dac_sym(1:Ldac:end) = IF_sym;

%% DMA
dma_sym = filter(DMA_new,dac_sym);
group_delay = 11;
dma_sym = dma_sym(group_delay+1:end);

%% DMA (Receiver)
dma_received = filter(DMA_new,dma_sym);
dma_received = dma_received(group_delay+1:end);

%%  ADC
Ladc = Ldac;
adc_received = dma_received(1:Ladc:end);

%% IF demodulation
m = 0:length(adc_received)-1;
IF_received = adc_received .* exp(-1i*2*pi*f_IF/fs*m);

%% Pulse shaping
srrc_received = conv(IF_received,SRRC_filter(trun,M,alpha),'same');

%% Downsampling
Ld = Lu;
downsample_received = srrc_received(1:Ld:end);

%% Phase
phase_received = unwrap(angle(downsample_received));
phase_diff = zeros(1,length(phase_received));
phase_diff(1) = phase_received(1);
for m=2:length(phase_received)
    phase_diff(m) = phase_received(m)-phase_received(m-1);
end

%% Mapping
bit = [0, 0, 1, 1; 0, 1, 1, 0];
mapper = [pi/4; 3*pi/4; -3*pi/4; -pi/4];
map_dist = abs(phase_diff-mapper);
[mini,index] = min(map_dist);
recover_bit = bit(:,index);
recover = reshape(recover_bit,1,length(recover_bit)*2);

%% Spectrum mask
[px,f] = pwelch(dma_sym,[],[],[],1);

figure(1);
p1 = (f_IF-2.5*1e6)/(1*1e6*Lu*Ldac);
p2 = (f_IF-1.5*1e6)/(1*1e6*Lu*Ldac);
p3 = (f_IF-1*1e6)/(1*1e6*Lu*Ldac);
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

%% EVM
symbol_constellation = exp(1i*symbol);
phase_diff_con = exp(1i*phase_diff);
EVM = 10*log10(mean(abs(symbol_constellation-phase_diff_con).^2) / mean(abs(symbol_constellation).^2));
disp(['The EVM for a noiseless channel is ',num2str(EVM),' dB']);
