f=1e7;    % frequency
d_1=500;    % distance
d_2=1000;    % 2nd path distance
g = 0.5;    %fading coeffecient
T=1/f/32; %sample rate
n=100;
t=0:T:n*T;
c=3*10^8; % light speed
l=c/f;    % wave length
mag1=zeros(1,length(.1:1:d_2));  % creat  1:d matrix mag1
mag2=zeros(1,length(.1:1:d_2));  % creat  1:d matrix mag2
k = 1;  % set counting variable
for r=1:1:d_2
      if r<=d_1 %first nest for r<=d_1
      c1=cos(2*pi*f*(t-r/c)) - cos(2*pi*f*(t-(2*d_1-r)/c));                 % no distance attenation for reflection wall signal
      c2=cos(2*pi*f*(t-r/c))/r - cos(2*pi*f*(t-(2*d_1-r)/c))/(2*d_1-r) - g*cos(2*pi*f*(t-(2*d_2-r)/c))/(2*d_2-r);   % with distance attenation for reflection wall signal
    
      mag1(k) = max(c1); % pick max value in a period T == magnitude
      mag2(k) = max(c2); % pick max value in a period T == magnitude 
      else       % d_1<r<d_2
      c1=cos(2*pi*f*(t-r/c)) - cos(2*pi*f*(t-(2*d_1-r)/c));                 % no distance attenation for reflection wall signal
      c2=cos(2*pi*f*(t-r/c))/r - g*cos(2*pi*f*(t-(2*d_2-r)/c))/(2*d_2-r);   % with distance attenation for reflection wall signal
    
      mag1(k) = max(c1); % pick max value in a period T == magnitude
      mag2(k) = max(c2); % pick max value in a period T == magnitude 
      end  
      k=k+1; 
      
end
figure(1)
plot(d_2/4:d_2,mag1(d_2/4:end));grid;
figure(2)
plot(d_2/4:d_2,mag2(d_2/4:end))