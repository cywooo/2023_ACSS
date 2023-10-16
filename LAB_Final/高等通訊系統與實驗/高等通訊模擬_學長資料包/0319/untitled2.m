f=1e7;    % frequency
d=200;    % distance
T=1/f/32; %sample rate
n=100;
t=0:T:n*T;
c=3*10^8; % light speed
l=c/f;    % wave length
mag1=zeros(1,length(.1:1:d));  % creat  1:d matrix mag1
mag2=zeros(1,length(.1:1:d));  % creat  1:d matrix mag2
k = 1;  % set counting variable
for r=1:1:d
    c1=cos(2*pi*f*(t-r/c)) - cos(2*pi*f*(t-(2*d-r)/c));             % no distance attenation for reflection wall signal
    c2=cos(2*pi*f*(t-r/c))/r - cos(2*pi*f*(t-(2*d-r)/c))/(2*d-r);   % with distance attenation for reflection wall signal
    
    mag1(k) = max(c1); % pick max value in a period T == magnitude
    mag2(k) = max(c2); % pick max value in a period T == magnitude 
    k=k+1; 
end 
figure(1)
plot(d/4:d,mag1(d/4:end));grid;
figure(2)
plot(d/4:d,mag2(d/4:end))