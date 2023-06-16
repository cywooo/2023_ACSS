sr = 10^6;
% ADC = 16*10^6;
% DMA = 64*10^6;
uf1 = 16;           % uf1 =   ADC(16M)/SR(1M) 
uf2 = 4;            % uf2 =   DMA(64M)/ADC(16M)

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

ffupup = filter(Hd,fupup);                   % IIR

figure(1)
subplot(5,1,1);stem(abs(fft(up)));title('up sampling in freq domain')
subplot(5,1,2);stem(abs(fft(fup)));title('SRRC Filter (Digital pulse shaping)')
subplot(5,1,3);stem(abs(fft(fupup)));title('2nd up sampling')
subplot(5,1,4);stem(abs(fft(ffupup)));title('After IIR Filter')
subplot(5,1,5);stem(ffupup);title('Origin tx signal in time domain')

%%
NPW=0.001;
n = sqrt(NPW)*randn(1,length(ffupup));
r = ffupup + n;

%% up conversion
afc = 16*10^6;   %analog carrier freq   (16M)
fc = afc/sr /(uf1*uf2);

x = real(r.*exp(1i*2*pi*(fc)*(1:size(ffupup,2))));   % x(m)

%% down conversion

xx = filter(IFF,x);

xe = xx.*cos((2*pi*(3/16)*(1:size(xx,2))) );

xef = filter(Hd,xe);                  % DMA filter   IIR
adjust =0;
xefd = xef(1+adjust:uf2:end);                %after DMA down  

XX = xefd.*exp( -1*1i *(2*pi*(1/4)*(1:size(xefd,2))));

xefdf = conv(XX,h,'same');          % Digital filter    SRRC
xefdfd = xefdf(1:uf1:end);            % = b(n)



figure(2)

subplot(5,1,1);stem(fftshift(abs(fft(x))));
subplot(5,1,2);stem(fftshift(abs(fft(xx))));title('先濾中頻')
subplot(5,1,3);stem(fftshift(abs(fft(xe))));title('先調到中頻')
subplot(5,1,4);stem(fftshift(abs(fft(xef))));title('在中頻濾掉')
subplot(5,1,5);stem(fftshift(abs(fft(xefdf))));


delay =1;
figure(3)
%stem(fftshift(abs(fft(cos(2*pi*(fc)*(1:size(ffupup,2)))))));
stem([BPSK(1:end-delay)' 23*xefdfd(1+delay:end)']);
% stem(fftshift(abs(fft(r))));

