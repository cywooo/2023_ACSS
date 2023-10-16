
n=10000;
m=10000;
x=sum(rand(m,n)-0.5);
y=(mean(x.^4)/(mean(x.^2))^2)-3;
hist(x);
