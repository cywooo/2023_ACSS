close all;
sr = 10^6;
uf1 = 16;
uf2 = 2;

%% BPSK
sn = 200;
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
r = ffupup;

%% up conversion
afc = 8*10^6;   %analog carrier freq
fc = afc/sr /(uf1*uf2);    % 1/4
g = 1.2;
phase = pi/6;
xi = sqrt(2)*real(r).*cos(2*pi*(fc)*(1:size(r,2)) );   % x(m)
xq = -g*sqrt(2)*imag(r).*sin((2*pi*(fc)*(1:size(r,2)) + phase) );

x = xi + xq;

%% down conversion

yi = sqrt(2)*x.*cos(2*pi*(fc)*(1:size(r,2)) );
yq = (-1)*sqrt(2)*x.*sin(2*pi*(fc)*(1:size(r,2)) );

yif = filter(LPF,yi); 
yqf = filter(LPF,yq); 

y = yif + 1i*yqf;
yf = filter(Hd,y);

delay =0;
yfd = yf(1+delay:uf2:end); 

yfdf = conv(yfd,h,'same'); 

yfdfd = yfdf(1:uf1:end);  

figure(2)
subplot(4,1,1);stem(y);title('Y_I freq')
subplot(4,1,2);stem(yf);title('Filter Y_I freq')
subplot(4,1,3);stem(yfd);title('Y_Q freq r')
subplot(4,1,4);stem(yi);title('Filter Y_Q freq')

% delay=2;
%figure(3)
scatterplot(yfdfd);



