sr = 10^6;
uf1 = 16;           % uf1 =   ADC(16M)/SR(1M) 
uf2 = 4;            % uf2 =   DMA(64M)/ADC(16M)
phase = pi/6;
g = 1;
%% QPSK
sn = 2000;
BPSK = randi([0,1],2,sn);
BPSK(BPSK==0)=-1;
QPSK = BPSK(1,:) + 1j*BPSK(2,:);
%SRRC
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

xi = sqrt(2)*real(r).*cos(2*pi*(fc)*(1:size(r,2)) );   % x(m)
xq = -g*sqrt(2)*imag(r).*sin((2*pi*(fc)*(1:size(r,2)) + phase) );



% x = xi + xq;

% x = real(r.*exp(1i*2*pi*(fc)*(1:size(ffupup,2))));   % x(m)

% xlp = filter(LPF,x);
%% down conversion
% yi = sqrt(2)*x.*cos(2*pi*(fc)*(1:size(r,2)) );
% yq = (-1)*sqrt(2)*x.*sin(2*pi*(fc)*(1:size(r,2)) );
% yif = filter(LPF,yi); 
% yqf = filter(LPF,yq); 
% y = yif + 1i*yqf;

xii = filter(IFF,xi); 
xqq = filter(IFF,xq); 

%xx = xlp;
% xx = filter(IFF,x);             % IF 
%%
%xe = xx.*cos((2*pi*(3/16)*(1:size(xx,2))) );

yi = sqrt(2)*xii.*cos(2*pi*(3/16)*(1:size(xii,2)) );
yq = (-1)*sqrt(2)*xqq.*sin(2*pi*(3/16)*(1:size(xqq,2)) );

yif = filter(LPF,yi); 
yqf = filter(LPF,yq); 
xe = yif + 1i*yqf;

xef = filter(Hd,xe);                
xefd = xef(1:uf2:end);              
%%
XX = xefd.*exp( -1*1i *(2*pi*(1/4)*(1:size(xefd,2))));
%%

xefdf = conv(XX,h,'same');          % Digital filter    SRRC
xefdfd = xefdf(1:uf1:end);            % = b(n)

scatterplot(xefdfd);

