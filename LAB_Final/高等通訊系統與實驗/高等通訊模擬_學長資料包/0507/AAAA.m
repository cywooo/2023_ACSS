sr = 10^6;
uf1 = 16;
uf2 = 2;

%% BPSK
sn = 20;
BPSK = randi([0,1],1,sn);
BPSK(BPSK==0)=-1;
%SRRC
%% digital pulse shaping 
h = rcosdesign(0.25,10,16,'sqrt');
up = zeros(1,uf1*sn);
up(1:uf1:end) = BPSK;                        % up 1
fup = conv(up,h,'same');                     % SRRC

fupup = zeros(1,uf2*size(fup,2));
fupup(1:uf2:end) = fup;                      % up 2
%ffupup = conv(fupup,h,'same');
ffupup = filter(Hd,fupup);                   % IIR
%ffupup = filter(SOS(1,1:3),SOS(1,5:6),fupup);

figure(1)
subplot(5,1,1);stem(ffupup);title('Origin signal sequence spectrum')
subplot(5,1,2);stem(BPSK);title('up sampling')
subplot(5,1,3);stem(abs(fft(fup)));title('SRRC Filter (Digital pulse shaping)')
subplot(5,1,4);stem(abs(fft(fupup)));title('2nd up sampling')
subplot(5,1,5);stem(abs(fft(ffupup)));title('After IIR Filter')

%%
r = ffupup;

%% up conversion
afc = 8*10^6;   %analog carrier freq
fc = afc/sr /(uf1*uf2);

x = real(r.*exp(1i*2*pi*(fc)*(1:size(ffupup,2))));   % x(m)

%% down conversion
xe = x.*exp(-1*1i*2*pi*(fc)*(1:size(ffupup,2)));

%xef = filter(SOS(1,1:3),SOS(1,5:6),xe);
xef = filter(Hd,xe);                  % DMA filter   IIR
xefd = xef(1:uf2*uf1:end);                %after DMA down  


figure(2)
subplot(4,1,1);stem(fftshift(abs(fft(x))));
subplot(4,1,2);stem(fftshift(abs(fft(xe))));
%subplot(5,1,2);stem(fftshift(abs(fft(cos(2*pi*(fc)*(1:size(ffupup,2)))))));
subplot(4,1,3);stem(fftshift(abs(fft(xef))));
subplot(4,1,4);stem(xef);
% subplot(5,1,5);stem(fftshift(abs(fft(xefdfd))));
delay=2;
figure(3)
stem([BPSK(1:end-delay)' 15*xefd(1+delay:end)'])

