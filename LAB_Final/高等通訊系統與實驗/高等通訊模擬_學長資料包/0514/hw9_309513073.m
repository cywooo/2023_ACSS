sr = 10^6;
uf1 = 16;     
uf2 = 4;    
CFO = 0;
g = 0;
phase = 0;
%% BPSK
sn = 20;
BPSK = randi([0,1],2,sn);
BPSK(BPSK==0)=-1;
QPSK = (BPSK(1,:) + 1j*BPSK(2,:) ) ;
%% digital pulse shaping 
h = rcosdesign(0.25,10,16,'sqrt');
up = zeros(1,uf1*sn);
up(1:uf1:end) = QPSK;                        % up 1
fup = conv(up,h,'same');                     % SRRC
fupup = zeros(1,uf2*size(fup,2));
fupup(1:uf2:end) = fup;                      % up 2
ffupup = filter(Hd,fupup);                   % IIR

%%
NPW=0;
n = sqrt(NPW)*randn(1,length(ffupup));
r = ffupup + n;

%% up conversion
afc = 16*10^6;   %analog carrier freq   (16M)
fc = afc/sr /(uf1*uf2);

xi = sqrt(2)*real(r).*cos(2*pi*(fc)*(1:size(r,2)) );  
xq = -g*sqrt(2)*imag(r).*sin((2*pi*(fc)*(1:size(r,2)) + phase) );
x = xi + xq;
%% down conversion

xx = filter(IFF,x);
xe = xx.*cos((2*pi*(3/16)*(1:size(xx,2))) );
xef = filter(Hd,xe);                  % DMA filter   IIR
adjust =0;
xefd = xef(1+adjust:uf2:end);                %after DMA down  
XX = xefd.*exp( -1*1i *(2*pi*(1/4)*(1:size(xefd,2))));
xefdf = conv(XX,h,'same');          % Digital filter    SRRC
xefdfd = xefdf(1:uf1:end);            % = b(n)


delay =1;
figure(3)
%stem(fftshift(abs(fft(cos(2*pi*(fc)*(1:size(ffupup,2)))))));
stem([BPSK(1:end-delay)' 23*xefdfd(1+delay:end)']);
stem(fftshift(abs(fft(r))));
% scatterplot(yfdfd);
