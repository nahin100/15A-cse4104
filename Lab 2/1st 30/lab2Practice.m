t = 0:0.01:2; % sample points from 0 to 2 in steps of 0.01
x = sin(2*pi*t); % Evaluate sin(2 pi t)
plot(t,x,'b'); % Create plot with blue line
xlabel('t in sec'); ylabel('x(t)'); % Label axis
title('Plot  of sin(2\pi t)'); % Title plot

n = 0:1:40; % sample index from 0 to 20
x = sin(0.1*pi*n); % Evaluate sin(0.2 pi n)
Hs = stem(n,x,'b','filled'); % Stem-plot with handle Hs
set(Hs,'markersize',4); % Change circle size
xlabel('n'); ylabel('x(n)'); % Label axis
title('Stem Plot of sin(0.2 \pi n)'); % Title plot

[x,n] = stepseq(0,-5,10);

Hs = stem(n,x,'b','filled'); % Stem-plot with handle Hs
xlabel('n'); ylabel('x(n)'); % Label axis
title('Stem Plot of sin(0.2 \pi n)');

set(Hs,'markersize',4); % Change circle size
set(gca, 'ylim', [-2 +2]);


n = [0:100]; x = (1.2).^n;
Hs = stem(n,x,'b','filled');

r=0.9;
theta=pi/10;
n = [0:100]; rn = r.^n; phin=mod(theta.*n, 2*pi) 
figure(1)
subplot(2,1,1), stem(n,rn,'b','filled');
subplot(2,1,2), stem(n,phin,'b','filled'), set(gca, 'ylim', [0 +2*pi]);

r=0.9;
theta=pi/10;
n = [0:100]; 
XR = (r.^n).*cos(theta.*n); 
XI = (r.^n).*sin(theta.*n);
figure(1)
subplot(2,1,1), stem(n,XR,'b','filled');
subplot(2,1,2), stem(n,XI,'b','filled');

n = [0:100]; 
x = 3*cos((pi/16)*n+(pi/13));
stem(n,x,'b','filled'); % Stem-plot with handle Hs

n=[-1,0,1,2,3];
x=[5,4,3,2,1];

m = -fliplr(n);
m1 = min([m,n]); m2 = max([m,n]); m = m1:m2;
nm = n(1)-m(1); n1 = 1:length(n);

x1 = [zeros(1,nm) x];
%x1(n1+nm) = x; 
x = x1;

xe = 0.5*(x + fliplr(x)); xo = 0.5*(x - fliplr(x));

figure(1)
subplot(2,1,1), stem(m,xe,'b','filled');
subplot(2,1,2), stem(m,xo,'b','filled');


n = [-10:9]; x = [5,4,3,2,1];
xtilde = x' * ones(1,4);
xtilde = (xtilde(:))';
subplot(2,2,4); stem(n,xtilde); title('Sequence in Problem 2.1d')
xlabel('n'); ylabel('xtilde(n)');

subplot(2,1,1); % Two rows, one column, first plot
plot(t,x,'b'); % Create plot with blue line

subplot(2,1,2); % Two rows, one column, second plot
Hs = stem(n,x,'b','filled'); % Stem-plot with handle Hs

%Signal shift
n1=input('Enter the amount to be delayed');
n2=input('Enter the amount to be advanced');
n=-2:2;
x=[-2 3 0 1 5];
subplot(3,1,1);
stem(n,x);
title('Signal x(n)');
m=n+n1;
y=x;
subplot(3,1,2);
stem(m,y);
title('Delayed signal x(n-n1)');
t=n-n2;
z=x;
subplot(3,1,3);
stem(t,z);
title('Advanced signal x(n+n2)');

%Signal
theta = 0:pi/100:2*pi;
y = sin(10*theta);
z = sin(10*theta)+cos(10*theta);

figure(1);
subplot(2,1,1), plot(x,y);
subplot(2,1,2), plot(x,z);


title('Two sine plot'),
xlabel('x data'),
ylabel('Amplitude');

set(gca, 'xlim', [0 2*pi]);
set(gca, 'ylim', [-3 +3]);

set(gca, 'xtick', [0.0 0.25 0.5 0.75 1.0]*2*pi);
set(gca, 'xticklabel', [{'0'} {'\pi/2'} {'\pi'} {'3\pi/2'} {'2\pi'}])

grid('on');
grid('minor');


%Practice 2nd Lab
t = 0:0.01:2;
x = sin(2*pi*t);
plot(t,x,'b');
xlabel('t is sec');
ylabel('x(t)');
title('plot of sin(2\pit');

set(gca, 'xlim', [0 2]);
set(gca, 'ylim', [-2 +2]);

n=0:1:40;
x=sin(0.2*pi*n);
Hs = stem(n,x,'b', 'filled');
set(Hs, 'markersize', 4);
title('stemplot of sin(2\pin)');

n1=-5;
n2=5;

n=n1:n2;
x=zeros(1,length(n));

for i=1:length(n)
    if n(i)==0
        x(i)=1
    end
end

stem(n,x);

%discrete sinusoidal signal for different frequency

n=0:1:100;
x=sin(2*pi*0.01*n);
subplot(4,1,1), stem(n,x,'b', 'filled');

n=0:1:100;
x=sin(2*pi*0.03*n);
subplot(4,1,2), stem(n,x,'b', 'filled');

n=0:1:100;
x=sin(2*pi*0.05*n);
subplot(4,1,3), stem(n,x,'b', 'filled');


n=[0:100];
r=0.9;
theta=pi/100;

xR=(r.^n).*cos(theta.*n);
xI=(r.^n).*sin(theta.*n);

figure(1);
subplot(2,1,1), stem(n,xR);

n=[-10:9];
x=[5,4,3,2,1];

xtilde=x'*ones(1,4);
xtilde=xtilde(:)';
stem(n,xtilde);

n=[-1,0,1,2,3];
x=[5,4,3,2,1];

m = -fliplr(n);
y=fliplr(x);

subplot(2,1,1),stem(n,x);
subplot(2,1,2),stem(m,y);

n=[-1,0,1,2,3];
x=[5,4,3,2,1];

m=-fliplr(n);
m1=min([n,m]); 
m2 = max([m,n]);

m=m1:m2;

nm=n(1)-m(1);
x1 = [zeros(1,nm), x];

x=x1;
xe=0.5*(x+fliplr(x));
xo=0.5*(x-fliplr(x));

subplot(3,1,1), stem(m,x,'b', 'filled');

subplot(3,1,2), stem(m,xe,'b', 'filled');

subplot(3,1,3), stem(m,xo,'b', 'filled');

x=[-2 3 0 1 5];
n=-2:2;

n1=2;

m=n+n1;
figure(1);
stem(n,x);
figure(2);
stem(m,x);

x1=[5,4,6,7,8];
n1=[-3,-2,-1,0,1];

x2=[10,11,12,13,14];
n2=[0,1,2,3,4];

m1 = min([n1,n2]);
m2 = max([n1,n2]);

m=m1:m2;

for i=m1:m2
  
end








