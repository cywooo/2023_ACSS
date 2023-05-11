clear all; clc; close all;

symbol_rate = 1e6;
DAC_rate = 16e6;
ANL_rate = 64e6;
IF_freq = 2e6;


M_1 = DAC_rate/symbol_rate;
M_2 = ANL_rate/DAC_rate;

SNR = [-10:5:40];
BT_ = [0.4 0.5 0.6]; %[0.3 0.4 0.5 0.6 0.7 0.8];
fd_ = [1.3 1.5 1.7]* symbol_rate;
BER = zeros(3,length(SNR));
%% HW
for d = [1:3]%[1:6] % 
    for g = [1:length(SNR)]
        signal_length = 1e3 +20;
        % signal generate------
        s_1 = [zeros(1,10) sign(randn(1,signal_length-20)) zeros(1,10)];
        signal_power = sum(abs(s_1(11:end-10)).^2)/(signal_length-20);
        % up sample by 16-------
        s_01_up_1 = up_sample(M_1,s_1);
        % Gussian filter ---------
        
        fd = 0.15 * symbol_rate; %fd_(d) %0.15 * symbol_rate; %150e3; %
        BT = BT_(d); %0.5; %
        B = BT*fd;

        GC_1 = 1;
        Gn_1 = [-floor(M_1*2-1):floor(M_1*2-1)];
        Gussian_filter = GC_1*exp(-1*((2*pi^2)/log(2))*(BT./ M_1)^2.*Gn_1.^2);
        
        if g==1
            figure;
            stem(Gn_1,Gussian_filter);
            title_text = "Gussian filter, n domain, M = "+num2str(M_1)+" ,BT = "+num2str(BT);
            title(title_text,"fontsize",12);
        end
        
        s_02_Gussianf = conv(s_01_up_1,Gussian_filter);
        s_02_Gussianf = s_02_Gussianf([floor((length(s_02_Gussianf)-length(s_01_up_1))/2)+1 :...
                    floor((length(s_02_Gussianf)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning
        % samation----------------
        %sum_s_02_Gussianf = s_02_Gussianf * triu(ones(length(s_02_Gussianf)));
        for k = [1:length(s_02_Gussianf)]
            sum_s_02_Gussianf(k) = sum(s_02_Gussianf([1:k]));
        end

        s_03_samation = exp(1i*2*pi*fd*(1/symbol_rate)*sum_s_02_Gussianf);

        % load IF frequency-------
        s_04_IF = s_03_samation .* exp(1i*2*pi*IF_freq/DAC_rate.*[1:length(s_03_samation)]);
        s_04_IF_real = real(s_04_IF);

        % up sample to analog-------
        s_05_up_analog = up_sample(M_2,s_04_IF_real);  

        s_05_up_analog_FeqDomain = fftshift(fft(s_05_up_analog));

        % LPF-----------------------
        s_06_toChannel = filter(Lab12_demo_IIR,s_05_up_analog);

        s_06_toChannel_FeqDomain = fftshift(fft(s_06_toChannel));
        

        % Channel-------------------
        noise_power = signal_power/(10^(SNR(g)/10));
        s_07_Channel = s_06_toChannel+randn(1,length(s_06_toChannel))*sqrt(noise_power);
        % LPF-----------------------
        s_08_receivedLPF = filter(Lab12_demo_IIR,s_07_Channel);
        % down sample--------------
        delay = 2;
        s_09_down = s_08_receivedLPF([delay:M_2:length(s_08_receivedLPF)]);

        % demodulate IF frequency-------
        s_10_deIFmodu = s_09_down .* exp(-1i*2*pi*IF_freq/DAC_rate.*[1:length(s_09_down)]);

        % SRRC----------------------
        %=========
        M_srrc = M_1;
        n_srrc = [-2*M_srrc:2*M_srrc] + 1e-6;% Avoid Singularity
        a = 0.5 ;% Avoid Singularity

        A_ = cos((1+a).*pi.*n_srrc./M_srrc);
        B_ = M_srrc.*sin((1-a).*pi.*n_srrc./M_srrc)./(4.*a.*n_srrc);
        C_ = 1-(4.*a.*n_srrc./M_srrc).^2;
        SRRC_n = (4.*a./pi).*(A_+B_)./C_;
        %=========
        s_11_SRRC = conv(s_10_deIFmodu,SRRC_n);
        s_11_SRRC = s_11_SRRC([floor((length(s_11_SRRC)-length(s_10_deIFmodu))/2)+1 :...
                    floor((length(s_11_SRRC)-length(s_10_deIFmodu))/2)+length(s_10_deIFmodu)]);% aligning
        % take phase------------------
        s_12_phase = unwrap(angle(s_11_SRRC))/2/pi/fd*DAC_rate/M_1;
        
        if (g==1 && d==1)
        figure;
        plot([1:length(s_12_phase)],s_12_phase);
        title_text = "phase";
        title(title_text,"fontsize",12);
        end 
        
        % diff-----------------------
        s_13_Diff = s_12_phase - [0 s_12_phase(1:end-1)] ;

        % Gussian filter ---------
        s_14_Gussianf = conv(s_13_Diff,Gussian_filter);
        s_14_Gussianf = s_14_Gussianf([floor((length(s_14_Gussianf)-length(s_13_Diff))/2)+1 :...
                    floor((length(s_14_Gussianf)-length(s_13_Diff))/2)+length(s_13_Diff)]);% aligning

        % down sample by 16-------------
        delay = 11;
        s_15_reDig = s_14_Gussianf([delay:M_1:length(s_14_Gussianf)]);

        s_16_compare = sign(s_15_reDig([11:11+signal_length-21]));
        s_1_compare = s_1([11:11+signal_length-21]);

        BER(d,g) = sum(abs(s_16_compare-s_1_compare))/2/(signal_length-20);

    end
end
%{
    BT_ = [0.4 0.5 0.6];
    fd_ = [1.3 1.5 1.7]* symbol_rate;
%}



figure;
plot(SNR,BER(1,:));
hold on;
plot(SNR,BER(2,:));
plot(SNR,BER(3,:));
%plot(SNR,BER(4,:));
%plot(SNR,BER(5,:));
%plot(SNR,BER(6,:));
hold off;
legend("BT=0.4","BT=0.5","BT=0.6");
%legend("BT=0.3","BT=0.4","BT=0.5","BT=0.6","BT=0.7","BT=0.8");
title_text = "Bit Error Rate in different BT, HW";

%legend("f_d=130kHZ","f_d=150kHZ","f_d=170kHZ");
%title_text = "Bit Error Rate in different f_d, HW";
xlabel("SNR");
ylabel("BER");
title(title_text,"fontsize",12);

