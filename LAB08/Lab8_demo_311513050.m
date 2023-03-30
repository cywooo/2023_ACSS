clear all; clc; close all;

M_1 = 4;
M_2 = 4;

n_1 = [-3*M_1:3*M_1] + 1e-6;% Avoid Singularity
n_2 = [-3*M_2:3*M_2] + 1e-6;% Avoid Singularity

%% Pulse fomula SRRC
a = 0.5+1e-6 ;% Avoid Singularity

% n domain SRRC
SRRC_1n = (4.*a./pi).*(cos((1+a).*pi.*n_1./M_1)+ M_1.*sin((1-a).*pi.*n_1./M_1)./(4.*a.*n_1))...
            ./(1-(4.*a.*n_1./M_1).^2);
        
SRRC_2n = (4.*a./pi).*(cos((1+a).*pi.*n_2./M_2)+ M_2.*sin((1-a).*pi.*n_2./M_2)./(4.*a.*n_2))...
            ./(1-(4.*a.*n_2./M_2).^2);
        
%% practice 1
% 1.Conduct SRRC pulse shaping for a QPSK sequence (the upsampling factor is 64).
qpsk_length = 1e2;
s_qPSK = sign(randn(1,qpsk_length));


% upsampling factor is 4
s_up_1 = up_sample(M_1,s_qPSK); % rate = symbol_rate*M
% SRRC
s_SRRC_1 = conv(s_up_1,SRRC_1n);
s_SRRC_1 = s_SRRC_1([floor((length(s_SRRC_1)-length(s_up_1))/2)+1 :...
            floor((length(s_SRRC_1)-length(s_up_1))/2)+length(s_up_1)]);% aligning

s_up_2 = up_sample(M_2,s_SRRC_1);       
s_SRRC_2 = conv(s_up_2,SRRC_2n);
s_SRRC_2 = s_SRRC_2([floor((length(s_SRRC_2)-length(s_up_2))/2)+1 :...
            floor((length(s_SRRC_2)-length(s_up_2))/2)+length(s_up_2)]);% aligning
        
s_out = conv(s_SRRC_2,SRRC_2n);
s_out = s_out([floor((length(s_out)-length(s_SRRC_2))/2)+1 :...
            floor((length(s_out)-length(s_SRRC_2))/2)+length(s_SRRC_2)]);% aligning
        
%% practice 2

%% practice 3

% show SRRC shape
figure;
plot(n/M,SRRC_n);
xlabel("n/M");
title_text = "SRRC shape, a ="+ num2str(a);
title(title_text,"fontsize",12);

% show s QPSK
figure;
stem([1:length(s_qPSK)],real(s_qPSK),"O");
hold on;
stem([1:length(s_qPSK)],imag(s_qPSK),"+");
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
