
%%

sr = 10^6;
uf1 = 16;
uf2 = 2;

%% BPSK
% sn = 1;
% BPSK = randi([0,1],1,sn);
% BPSK(BPSK==0)=-1;
N=32;
Impluse = [1,zeros(1,N-1)];

%% digital pulse shaping 
h = rcosdesign(0.25,10,16,'sqrt');

up = zeros(1,uf1*N);
up(1:uf1:end) = Impluse;                        % up 1
fup = conv(up,h,'same');                     % SRRC

fupup = zeros(1,uf2*size(fup,2));
fupup(1:uf2:end) = fup;                      % up 2
ffupup = filter(Hd,fupup);                   % IIR   xb


%% up conversion
afc = 8*10^6;   %analog carrier freq     fc = 8M
fc = afc/sr /(uf1*uf2);                 %1/4

x = real(ffupup.*exp(1i*2*pi*(fc)*(1:size(ffupup,2))));   % x(m)
de = 31;
h = [1 zeros(1,de-1) 0.5];
y = conv(x,h);
y = y(1:end-de);

%% Rx
yc = y.*exp(-1*1i*2*pi*(fc)*(1:size(y,2)));

xef = filter(Hd,yc);                  % DMA filter   IIR
xefd = xef(1:uf2*uf1:end);  

delay=2;
figure(2)
stem([Impluse(1:end-delay)' 10*xefd(1+delay:end)'])










