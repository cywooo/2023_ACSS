clear all; clc; close all;

M = 32;
n = [-3*M:3*M] + 1e-6;
t = linspace(-3/2,3/2,1e4);
f = linspace(-3,3,length(t));

%% Patice01
a = 0.9 ;
a = max(abs(a),1e-6);% avoid nan value
W = 0.5;

T = 1/(2*W);
f1 = W*(1-a);

% t domain
RC_t = (sin(pi*2.*W.*t)./(pi*2.*W.*t))...
        .*(cos(2*pi.*a.*W.*t)./((1-(4.*a.*W.*t).^2)+1e-5.*(abs(1-(4.*a.*W.*t).^2)==0)));
    
RC_f_th = 1/(2*W).*(0 <= abs(f)).*(abs(f) < f1)+...
        (1/(4*W)*(1-sin((pi*(abs(f)-W))/(2*W-2*f1)))).*(f1 <= abs(f)).*(abs(f) < 2*W-f1)+...
        0.*(f1 <= abs(f));

SRRC_t = (4.*a./pi).*(cos((1+a).*pi.*t./T)+T.*sin((1-a).*pi.*t./T)./(4.*a.*t))./(1-(4.*a.*t./T).^2);


% n domain
RC_n = (sin(pi*2.*W.*n.*T./M)./(pi*2.*W.*n.*T./M))...
        .*(cos(2*pi.*a.*W.*n.*T./M)./((1-(4.*a.*W.*n.*T./M).^2)));
    
SRRC_n = (4.*a./pi).*(cos((1+a).*pi.*n./M)+ M.*sin((1-a).*pi.*n./M)./(4.*a.*n))...
            ./(1-(4.*a.*n./M).^2);
        
RC_w = fftshift(fft(RC_n));
SRRC_w = fftshift(fft(RC_n));


figure(1);
plot(n/M,RC_n);
title_text = "RC n domain, a ="+num2str(a);
title(title_text,"fontsize",12);

figure(2);
plot(n/M,SRRC_n);
title_text = "SRRC n domain, a ="+num2str(a);
title(title_text,"fontsize",12);

figure(3);
plot(linspace(-1/2,1/2,length(SRRC_w)),abs(SRRC_w));
title_text = "SRRC n domain, a ="+num2str(a);
title(title_text,"fontsize",12);

SRRC_conv = conv_(SRRC_n,SRRC_n);

figure(4);
plot(n/M,RC_n,'-R');
hold on;
plot(n/M,SRRC_conv((length(SRRC_conv)/2-length(RC_n)/2+1):(length(SRRC_conv)/2+length(RC_n)/2))/M,'--B');
hold off;
legend("RC","conv. SRRC","fontsize",12);

%% Patice02
bpsk_length = 100;
s_BPSK = sign(randn(1,bpsk_length));
s_up = up_sample(M,s_BPSK);

s_pulse_RC = conv(s_up,RC_n);
s_RC_down = down_sample(M,s_pulse_RC);

s_up_sq = conv(s_up,ones(1,M));
s_up_SRRC_1 = conv(s_up_sq,SRRC_n);
s_pulse_SRRC = conv(s_up_SRRC_1,SRRC_n);

s_SRRC_down = down_sample(M,s_pulse_SRRC([floor((length(s_pulse_SRRC)-length(s_up))/2):end]));
s_SRRC_dmean = sum(abs(s_SRRC_down))/bpsk_length;

% show pulse shaping RC
figure(5);
stem([1:length(s_up)]+floor(length(n)/2),s_up);
hold on;
plot([1:length(s_pulse_RC)],s_pulse_RC);
hold off;
legend("s BPSK up x32","RC pulsed");
title_text = "s pulse shaping, RC";%+num2str(a)
title(title_text,"fontsize",12);
% show down-sample to pulse shaping RC
figure(6);
stem([1:length(s_BPSK)]+floor((length(s_RC_down)-length(s_BPSK))/2),s_BPSK,"XR");
hold on;
stem([1:length(s_RC_down)],s_RC_down,"+B");
hold off;
legend("s BPSK","RC pulse downsample")
title_text = "downsample to RC pulse shaped s";
title(title_text,"fontsize",12);

% show pulse shaping SRRC
figure(7);
stem([1:length(s_up)]+floor((length(s_pulse_SRRC)-length(s_up))/2),s_up);
hold on;
plot([1:length(s_pulse_SRRC)],s_pulse_SRRC/s_SRRC_dmean);
hold off;
legend("s BPSK up x32","SRRC pulsed");
title_text = "s pulse shaping, SRRC";%+num2str(a)
title(title_text,"fontsize",12);
% show down-sample to pulse shaping SRRC
figure(8);
stem([1:length(s_BPSK)],s_BPSK,"XR");%+length(n)/2-0.5
hold on;
stem([1:length(s_SRRC_down)],s_SRRC_down/s_SRRC_dmean,"oB");
hold off;
legend("s BPSK","RC pulse downsample")
title_text = "downsample to SRRC pulse shaped s";
title(title_text,"fontsize",12);

%% Patice03
s_BPSK_IIR = filter(lab7_IIR_filter,s_up);
s_BPSK_IIR = filter(lab7_IIR_filter,s_BPSK_IIR);

s_BPSK_IIR_d = down_sample(M,s_BPSK_IIR([M:end]));
s_BPSK_IIR_dmean = sum(abs(s_BPSK_IIR_d))/bpsk_length;

% show IIR pulse shaping
figure(9);
stem([1:length(s_up)]+M,s_up);
hold on;
plot([1:length(s_BPSK_IIR)],s_BPSK_IIR/s_BPSK_IIR_dmean);
hold off;
legend("s BPSK up x32","IIR filter pulsed");
title_text = "s pulse shaping, IIR filter";
title(title_text,"fontsize",12);
% show down-sample to pulse shaping SRRC
figure(10);
stem([1:length(s_BPSK)],s_BPSK,"XR");%+length(n)/2-0.5
hold on;
stem([1:length(s_BPSK_IIR_d)],s_BPSK_IIR_d/s_BPSK_IIR_dmean,"oB");
hold off;
legend("s BPSK","IIR filter pulse downsample")
title_text = "downsample to IIR filter pulse shaped s";
title(title_text,"fontsize",12);



