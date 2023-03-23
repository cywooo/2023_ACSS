clear all; clc; close all;


t = linspace(-3/2,3/2,1e4);
f = linspace(-3,3,length(t));

%% Patice01
a = 1/2 ;
a = max(abs(a),1e-6);% avoid nan value
W = 1;
T = 1/(2*W);
f1 = W*(1-a);

RC_t = (sin(pi*2.*W.*t)./(pi*2.*W.*t))...
        .*(cos(2*pi.*a.*W.*t)./((1-(4.*a.*W.*t).^2)+1e-5.*(abs(1-(4.*a.*W.*t).^2)==0)));

RC_f_th = 1/(2*W).*(0 <= abs(f)).*(abs(f) < f1)+...
        (1/(4*W)*(1-sin((pi*(abs(f)-W))/(2*W-2*f1)))).*(f1 <= abs(f)).*(abs(f) < 2*W-f1)+...
        0.*(f1 <= abs(f));

RC_f = fftshift(fft(RC_t));

SRRC = (4.*a./pi).*(cos((1+a).*pi.*t./T)+T.*sin((1-a).*pi.*t./T)./(4.*a.*t))./(1-(4.*a.*t./T).^2);

figure;
plot(t/T,RC_t);
title_text = "Raised Cosine time domain, a ="+num2str(a);
title(title_text,"fontsize",12);
figure;
plot(f/W,2.*W.*RC_f_th);
title_text = "Raised Cosine frequency domain (theorical), a ="+num2str(a);
title(title_text,"fontsize",12);
figure;
plot(f/W,abs(RC_f));
title_text = "Raised Cosine frequency domain (FFT), a ="+num2str(a);
title(title_text,"fontsize",12);
figure;
plot(t/T,SRRC);
title_text = "Squared Root Raised Cosine time domain, a ="+num2str(a);
title(title_text,"fontsize",12);

SRRC_conv = conv_(SRRC,SRRC);

figure;
plot([1:length(SRRC_conv)],SRRC_conv);%t/T
title_text = "Conv. SRRC, a ="+num2str(a);
title(title_text,"fontsize",12);

figure(5);
plot(t/T,RC_t,'-R');
hold on;
figure(5);
plot(t/T,SRRC_conv((length(SRRC_conv)/2-length(RC_t)/2+1):(length(SRRC_conv)/2+length(RC_t)/2))/1665,'--B');
hold off;
legend("RC","conv. SRRC","fontsize",12);
%% Patice02
bpsk_length = 1e3;
s_BPSK = sign(randn(1,bpsk_length));



%% Patice03






