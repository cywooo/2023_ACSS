
%load SOS;
SOS = [1,1.92664029874089,1,1,-1.82884566246877,0.850945380024170];
sr = 10^6;      %  Digital freq = Symbol rate

uf1 = 4;
uf2 = 8;

%% BPSK
sn = 20;
BPSK = randi([0,1],1,sn);
BPSK(BPSK==0)=-1;

%% digital pulse shaping                          up sampling & Digital filter

h = rcosdesign(0.25,10,6,'sqrt');
up = zeros(1,uf1*sn);
up(1:uf1:end) = BPSK;
fup = conv(up,h,'same');          % SRRC


%up sampling total = 32            up sampling & DMA filter
fupup = zeros(1,uf2*size(fup,2));
fupup(1:uf2:end) = fup;
%ffupup = conv(fupup,h,'same');
ffupup = filter(SOS(1,1:3),SOS(1,5:6),fupup);

figure(1)
subplot(5,1,1);stem(abs(fft(BPSK)));title('Origin signal sequence spectrum')
subplot(5,1,2);stem(abs(fft(up)));title('up sampling')
subplot(5,1,3);stem(abs(fft(fup)));title('SRRC Filter (Digital pulse shaping)')
subplot(5,1,4);stem(abs(fft(fupup)));title('2nd up sampling')
subplot(5,1,5);stem(abs(fft(ffupup)));title('After IIR Filter')
r = ffupup;
%% Rx
rf = filter(SOS(1,1:3),SOS(1,5:6),r);
rfd = rf(1:uf2:end);
rfdf = conv(rfd,h,'same');
rfdfd = rfdf(1:uf1:end);
% rfdfd(rfdfd>0)=1;
% rfdfd(rfdfd<0)=-1;
figure(2)
stem([BPSK(1:end)' 0.99*rfdfd(1:end)']);legend('BPSK','recover')

%% up conversion
afc = 8*10^6;   %analog carrier freq

fc = afc/sr /(uf1*uf2);
x = real(r.*exp(1i*2*pi*(fc)*(1:size(ffupup,2))));


figure(3)

stem(abs(fft(x)));title('Up conversion signal sequence spectrum')

