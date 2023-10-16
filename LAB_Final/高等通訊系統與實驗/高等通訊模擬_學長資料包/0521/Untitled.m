sr = 10^6;
uf1 = 16;
uf2 = 2;

%% BPSK
sn = 1;
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

% ffupup = xb                                  

%% up conversion
afc = 8*10^6;   %analog carrier freq     fc = 8M
fc = afc/sr /(uf1*uf2);
x = real(ffupup.*exp(1i*2*pi*(fc)*(1:size(ffupup,2))));   % x(m)

%% channel
phase=0;

% t = [1:size(ffupup,2)];
% a_b0 = 1*exp(-1*1i*(2*pi*fc*t+phase));
% a_b1 = 0.5*exp(-1*1i*(2*pi*fc*t+phase));
H=[1 ,zeros(30), 0.5];    % a_b(0) a_b(1)
r = conv(ffupup,H,'same'); 

figure(3)
stem(r);
%% down conversion
%xe = r.*exp(-1*1i*2*pi*(fc)*(1:size(ffupup,2)));
 
xef = filter(Hd,r);                  % DMA filter   IIR
xefd = xef(1:uf2*uf1:end);                %after DMA down  

% delay=4;
% figure(2)
% stem([BPSK(1:end-delay)' 15*xefd(1+delay:end)'])
