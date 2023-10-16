sr = 10^6;
IF = 4*10^6;
uf1 = 16;
uf2 = 2;


%% BPSK
sn = 20;
BPSK = randi([0,1],1,sn);
BPSK(BPSK==0)=-1;



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
r = ffupup;

%% up conversion
afc = 8*10^6;   %analog carrier freq
fc = afc/sr /(uf1*uf2);

x = real(r.*exp(1i*2*pi*(fc)*(1:size(ffupup,2))));   % x(m)

%% down conversion
xc = x.*exp((-1i*2*pi*(1/8)*(1:size(x,2))));

xcf = filter(Hd,xc);                  % DMA filter   IIR

xcfd = xcf(1:uf2:end);                %after DMA down  

xcfde = xcfd.*exp(-1*1i*2*pi*(1/4)*(1:size(xcfd,2)) );

xcfdef = conv(xcfde,h,'same');          % Digital filter    SRRC

xcfedfd = xcfdef(1:uf1:end);            % = b(n)


figure(2)
subplot(4,1,1);stem(fftshift(abs(fft(x))));
subplot(4,1,2);stem(fftshift(abs(fft(xc))));
subplot(4,1,3);stem(fftshift(abs(fft(xcf))));
subplot(4,1,4);stem(xcf);title('recover signal in time domain')

delay=3;
figure(3)
stem([BPSK(1:end-delay)' 15*xcfedfd(1+delay:end)'])


