t=0:0.01:2;
x=sin(2*pi*t);
plot(t,x);

%stem plot use
n=0:1:40;
x=sin(0.2*pi*n);
Hs=stem(n,x,'r','filled');
set(Hs,'markersize',4);
xlabel('n'); ylabel('x(n)');
title('stemplot of sin(0.2\pin)');

% Elementary signals (step signal)
n1=-5;
n2=5;

n=n1:n2;
x=zeros(1,length(n));

for i=1:length(n)
    if n(i)==0
        x(i)=1;
    end
end

stem(n,x);

%discrete sinusoid signal for different frequencies
n=0:1:100;
a=1;

x=a*sin(2*pi*0.01*n);
subplot(3,1,1); stem(n,x,'r','filled');

x=a*sin(2*pi*0.04*n);
subplot(3,1,2); stem(n,x,'r','filled');

x=a*sin(2*pi*0.09*n);
subplot(3,1,3); stem(n,x,'r','filled');

%Complex exponential
n=0:100;
r=0.9;
theta=pi/100;

xR=(r.^n).*cos(theta.*n);
xI=(r.^n).*sin(theta.*n);

subplot(2,1,1); stem(n,xR);
subplot(2,1,2); stem(n,xI);

%Creating periodic signals
x=[1,2,3,4,5,6];
n=0:23;

xtilde=x'*ones(1,4);
xtilde=xtilde(:)';

stem(n,xtilde,'b', 'filled');

% Signal fold x(n)
x=[1,2,3,4,5,6];
n=[-2,-1,0,1,2,3];

for i=1:length(n)
    m(i)=-n(length(n)-i+1);
end

for i=1:length(n)
    x1(i)=x(length(n)-i+1);
end

subplot(2,1,1); stem(n,x);
subplot(2,1,2); stem(m,x1);

% Even and odd component of any signal
n=[-1,0,1,2,3];
x=[5,4,3,2,1];

for i=1:length(n)
    m(i)=-n(length(n)-i+1);
end

m1 = min([m,n]);
m2 = max([m,n]);

m = m1:m2;
extend=n(1)-m(1);


x1 = [zeros(1,extend), x];

for i=1:length(m)
    x2(i)=x1(length(m)-i+1);
end

xe = 0.5*(x1+x2);
xo = 0.5*(x1-x2);

subplot(3,1,1); stem(m,x1);
subplot(3,1,2); stem(m,xe);
subplot(3,1,3); stem(m,xo);

% Signal shift
% Signal delay
n=[-2,-1,0,1,2,3];
x=[6,5,4,3,2,1];

k=2;

n1=n+k;

subplot(2,1,1); stem(n,x);
subplot(2,1,2); stem(n1,x);
% adjust index in a way that shift delay
% can be observed


% Signal addition
% Signal multiplication














































