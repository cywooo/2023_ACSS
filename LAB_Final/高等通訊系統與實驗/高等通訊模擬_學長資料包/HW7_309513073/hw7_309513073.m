uf = 64;
sr = 10^6;      %  Digital freq
fc = 8*10^6;
a = 0.25 ;
Ns = 100;  % # of sequence

%% SRRC set
M = uf;
n = [-128:127]-0.0001;
SRRC = (4*a/pi)* ( cos((1+a)*pi*n./M ) + (M*sin((1-a)*pi*n./M)  ./ (4*a*n) )  ) ./ (1-(4*a*n./M).^2);
%% QPSK sequence
QPSK = randi([0 1],2,Ns); % random signal 2*Ns
QPSK(QPSK==0) = -1;
QPSK_sequence = (1/sqrt(2)) *(QPSK(1,:) + 1j*QPSK(2,:)); % symbol with power 1

% BPSK = randi([0 1],1,Ns); %for practical realize
% BPSK(BPSK==0) = -1;
%% upsampling
du = zeros(1,Ns*uf);
du(1:uf:end)=QPSK_sequence;

%% filter
DMA = ones(1,64);
%fd = filter(Hd,du);,
fd = conv(DMA,du);

ffd = conv(SRRC,fd);
Np = size(ffd,2);
%% *exp
e = exp(1i*2*pi*(1/8)*[1:Np]);
X =ffd.*e;
x =real(X);
%% Fiigure
% figure(1)
% stem(QPSK_sequence);
% title('BPSK sequence')
% figure(2)
% stem(du);
% title('upsampling')
% figure(3)
% plot(fd);
% title('Filter (hold)')
% figure(4)
% plot(ffd);
% title('Filter')
figure(5)
stem(abs(fft(x)));
title('Re part Freq Domain')
% figure(6)
% stem(abs(fft(cos(2*pi*0.125*[1:6718]))));
% title('time domain')
% figure(7)
% plot(abs(fft(SRRC)));