close all;
sr = 10^6;
uf1 = 16;
uf2 = 2;
g = 1;
phase = pi/6;
%% BPSK
sn = 2000;
BPSK = randi([0,1],2,sn);
BPSK(BPSK==0)=-1;
QPSK = BPSK(1,:) + 1j*BPSK(2,:);

%% digital pulse shaping 
h = rcosdesign(0.25,10,16,'sqrt');
up = zeros(1,uf1*sn);
up(1:uf1:end) = QPSK;                        % up 1
fup = conv(up,h,'same');                     % SRRC

fupup = zeros(1,uf2*size(fup,2));
fupup(1:uf2:end) = fup;                      % up 2
ffupup = filter(Hd,fupup);                   % IIR

%%

%% up conversion
afc = 8*10^6;   %analog carrier freq
fc = afc/sr /(uf1*uf2);    % 1/4
xi = sqrt(2)*real(ffupup).*cos(2*pi*(fc)*(1:size(ffupup,2)) );   % x(m)
xq = -g*sqrt(2)*imag(ffupup).*sin((2*pi*(fc)*(1:size(ffupup,2)) + phase) );
x = xi + xq;



NPW=0.01;
n = sqrt(NPW)*randn(1,length(x));
r = x + n;
%% IF
xx = filter(IFF,x);

%% down conversion
yi = sqrt(2)*r.*cos(2*pi*(fc)*(1:size(r,2)) );
yq = (-1)*sqrt(2)*r.*sin(2*pi*(fc)*(1:size(r,2)) );

yif = filter(LPF,yi); 
yqf = filter(LPF,yq); 

y = yif + 1i*yqf;
yf = filter(Hd,y);

delay =0;
yfd = yf(1+delay:uf2:end); 

yfdf = conv(yfd,h,'same'); 

yfdfd = yfdf(1:uf1:end);  


scatterplot(yfdfd);



