clear all
close all

BT_list = [0.1:0.1:1];
BER = zeros(1,length(BT_list));
iteration = 1000;

for s=1:length(BT_list)
    BER_temp = 0;
    for iter=1:iteration
              %% BPSK
        N = 64;
        symbol = randi([0 1],1,N);
        symbol(symbol==0) = -1;

              %% Upsample
        Lu = 16;
        upsample_sym = zeros(1,length(symbol)*Lu);
        upsample_sym(1:Lu:end) = symbol;

              %% Gaussian filter
        BT = BT_list(s);
        M = Lu;
        trun = 5;
        gaussian_sym = conv(upsample_sym,Gaussian_filter(BT,M,trun),'same');

              %% Summation to phase (Integration)
        sf = zeros(1,length(gaussian_sym));
        for m=1:length(gaussian_sym)
            sf(m) = sum(gaussian_sym(1:m));
        end
        fd = 150*1e3;
        fb = 1*1e6;
        Tb = 1/(fb*M);
        sum_sym = exp(1i*2*pi*fd*Tb.*sf);

              %% IF
        f_IF = 2*1e6;
        fs = fb*Lu;
        m = 0:length(sum_sym)-1;
        IF_sym = real(sum_sym .* exp(1i*2*pi*f_IF/fs*m));

              %% DAC
        Ldac = 4;
        dac_sym = zeros(1,length(IF_sym)*Ldac);
        dac_sym(1:Ldac:end) = IF_sym;

              %% DMA
        dma_sym = filter(DMA,dac_sym);
        group_delay = 3;
        dma_sym = dma_sym(group_delay+1:end);
%         alpha = 0.3;
%         dma_sym = conv(dac_sym,SRRC_filter(trun,4,alpha),'same');

              %% DMA (Receiver)
        dma_receive = filter(DMA,dma_sym);
        dma_receive = dma_receive(group_delay+1:end);
%         alpha = 0.3;
%         dma_receive = conv(dma_sym,SRRC_filter(trun,4,alpha),'same');

              %% ADC
        Ladc = Ldac;
        adc_receive = dma_receive(1:Ladc:end);

              %% IF
        m = 0:length(adc_receive)-1;
        IF_receive = adc_receive .* exp(-1i*2*pi*f_IF/fs*m);

              %% SRRC filter
        alpha = 0.3;
        srrc_receive = conv(IF_receive,SRRC_filter(trun,M,alpha));
        delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
        srrc_receive = srrc_receive(delay+1:end-delay);

              %% Get phase (Differentiation)
        phase_sum = phase(srrc_receive) / (2*pi*fd*Tb);
        diff = zeros(1,length(phase_sum));
        diff(1) = phase_sum(1);
        for m = 2:length(phase_sum)
            diff(m) = phase_sum(m)-phase_sum(m-1);
        end

              %% Gaussian filter
        gaussian_receive = conv(diff,Gaussian_filter(BT,M,trun),'same');

              %% Downsample
        Ld = Lu;
        downsample_r = gaussian_receive(1:Ld:end);

              %% Detection
        detection = (downsample_r>0) - (downsample_r<0);
        
              %% BER
        BER_temp = BER_temp + sum(detection~=symbol);
        
    end
    BER(s) = BER_temp/(iteration*length(symbol));
end

%% Plot
figure(1);
plot(BT_list,BER,'--o');
xlabel('Bandwidth (BT)');
ylabel('BER');
title('BER with different bandwidth');
