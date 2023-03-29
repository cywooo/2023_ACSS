clear all; clc; close all;

M = 64;
n = [-3*M:3*M] + 1e-6;% Avoid Singularity
symbol_rate = 1e6;
carrier_freq_analog = 8e6;

%% Pulse fomula SRRC
a = 0.5+1e-6 ;% Avoid Singularity

% n domain SRRC
SRRC_n = (4.*a./pi).*(cos((1+a).*pi.*n./M)+ M.*sin((1-a).*pi.*n./M)./(4.*a.*n))...
            ./(1-(4.*a.*n./M).^2);
        
%% HW
% 1.Conduct SRRC pulse shaping for a QPSK sequence (the upsampling factor is 64).
qpsk_length = 1e2;
s_QPSK = sqrt(2).*sign(randn(1,qpsk_length)) + sqrt(2).*sign(randn(1,qpsk_length))*i;


% upsampling factor is 64
s_up = up_sample(M,s_QPSK);% rate = symbol_rate*M

% Use the practical DAC.
s_up_sq = conv(s_up,ones(1,M));
s_up_sq = s_up_sq([floor((length(s_up_sq)-length(s_up))/2)+1 :...
            floor((length(s_up_sq)-length(s_up))/2)+length(s_up)]);% aligning 
        
s_up_SRRC_1 = conv(s_up_sq,SRRC_n);
s_up_SRRC_1 = s_up_SRRC_1([floor((length(s_up_SRRC_1)-length(s_up_sq))/2)+1 :...
            floor((length(s_up_SRRC_1)-length(s_up_sq))/2)+length(s_up_sq)]);% aligning
        
s_up_SRRC_2 = conv(s_up_SRRC_1,SRRC_n);
s_up_SRRC_2 = s_up_SRRC_2([floor((length(s_up_SRRC_2)-length(s_up_SRRC_1))/2)+1 :...
            floor((length(s_up_SRRC_2)-length(s_up_SRRC_1))/2)+length(s_up_SRRC_1)]);% aligning

% Let the symbol rate be 1MHz, the carrier frequency be 8MHz. 
% Conduct the up-conversion operation in the equivalent digital domain.
w_c = carrier_freq_analog/symbol_rate/M*2*pi;
carrier = exp(i*w_c.*[1:length(s_up_SRRC_2)]);

s_out = real(s_up_SRRC_2 .* carrier);

S_up_SRRC_2 = fftshift(fft(s_up_SRRC_2));
S_out = fftshift(fft(s_out));

% show SRRC shape
figure;
plot(n/M,SRRC_n);
xlabel("n/M");
title_text = "SRRC shape, a ="+ num2str(a);
title(title_text,"fontsize",12);

% show s QPSK
figure;
stem([1:length(s_QPSK)],real(s_QPSK),"O");
hold on;
stem([1:length(s_QPSK)],imag(s_QPSK),"+");
hold off;
legend("Real","Imag","fontsize",12);
title_text = "s QPSK";
title(title_text,"fontsize",12);

% show s before carried
figure;
plot([1:length(s_up_SRRC_2)],real(s_up_SRRC_2));
hold on;
plot([1:length(s_up_SRRC_2)],imag(s_up_SRRC_2));
hold off;
legend("Real","Imag","fontsize",12);
title_text = "s before carrier";
title(title_text,"fontsize",12);

% show carrier
figure;
stem([1:length(carrier)],real(carrier));
title_text = "carrier";
title(title_text,"fontsize",12);

% Observe the up-converted spectrum to see if your design is correct.
figure;
plot(linspace(-1/2*symbol_rate*M,1/2*symbol_rate*M,length(S_out)),abs(S_out));
title_text = "s_o_u_t Frequency Domain";
title(title_text,"fontsize",12);
