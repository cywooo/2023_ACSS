clear all; clc; close all;

M_1 = 4;
M_2 = 4;

n_1 = [-3*M_1:3*M_1] + 1e-6;% Avoid Singularity
n_2 = [-3*M_2:3*M_2] + 1e-6;% Avoid Singularity

%% Pulse fomula SRRC
a = 0.5+1e-6 ;% Avoid Singularity

% n domain SRRC
SRRC_1_n = (4.*a./pi).*(cos((1+a).*pi.*n_1./M_1)+ M_1.*sin((1-a).*pi.*n_1./M_1)./(4.*a.*n_1))...
            ./(1-(4.*a.*n_1./M_1).^2);
        
SRRC_2_n = (4.*a./pi).*(cos((1+a).*pi.*n_2./M_2)+ M_2.*sin((1-a).*pi.*n_2./M_2)./(4.*a.*n_2))...
            ./(1-(4.*a.*n_2./M_2).^2);
        
%% practice 1
% 1.Conduct SRRC pulse shaping for a QPSK sequence (the upsampling factor is 64).
qpsk_length = 1e2;
s_qPSK = sign(randn(1,qpsk_length));


% upsampling factor is 4
s_01_up_1 = up_sample(M_1,s_qPSK); % rate = symbol_rate*M

% SRRC
s_02_SRRC_1 = conv(s_01_up_1,SRRC_1_n);
s_02_SRRC_1 = s_02_SRRC_1([floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+1 :...
            floor((length(s_02_SRRC_1)-length(s_01_up_1))/2)+length(s_01_up_1)]);% aligning

s_03_up_2 = up_sample(M_2,s_02_SRRC_1);       

s_04_SRRC_2 = conv(s_03_up_2,SRRC_2_n);
s_04_SRRC_2 = s_04_SRRC_2([floor((length(s_04_SRRC_2)-length(s_03_up_2))/2)+1 :...
            floor((length(s_04_SRRC_2)-length(s_03_up_2))/2)+length(s_03_up_2)]);% aligning
%{       
s_05_out = conv(s_04_SRRC_2,SRRC_2_n);
s_05_out = s_05_out([floor((length(s_05_out)-length(s_04_SRRC_2))/2)+1 :...
            floor((length(s_05_out)-length(s_04_SRRC_2))/2)+length(s_04_SRRC_2)]);% aligning
%}
        
%% practice 2
s_05_add_noise = s_04_SRRC_2 + randn(1,length(s_04_SRRC_2)).*sqrt(1);

s_06_received_SRRC2 = conv(s_05_add_noise,SRRC_2_n);
s_06_received_SRRC2 = s_06_received_SRRC2([floor((length(s_06_received_SRRC2)-length(s_05_add_noise))/2)+1 :...
            floor((length(s_06_received_SRRC2)-length(s_05_add_noise))/2)+length(s_05_add_noise)]);% aligning

s_07_down_1 = down_sample(M_2,s_06_received_SRRC2);

s_08_received_SRRC1 = conv(s_07_down_1,SRRC_1_n);
s_08_received_SRRC1 = s_08_received_SRRC1([floor((length(s_08_received_SRRC1)-length(s_07_down_1))/2)+1 :...
            floor((length(s_08_received_SRRC1)-length(s_07_down_1))/2)+length(s_07_down_1)]);% aligning

a_hat_09 = down_sample(M_1,s_08_received_SRRC1);

a_hat_normalized = a_hat_09/(sqrt(mean(a_hat_09.^2)));

figure;
stem([1:length(a_hat_normalized)],real(a_hat_normalized));
hold on;
stem([1:length(s_qPSK)],real(s_qPSK));
hold off;

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
