%Signal Addition and Multiplication
n1=[-1,0,1];
x1=[5,4,3];

n2=[-2,-1,0,1,2,3];
x2=[1,2,3,4,7,8,9];

m1 = min([n1,n2]);
m2 = max([n1,n2]);

m = m1:m2;

newx1= zeros(1, length(m));
newx2= zeros(1, length(m));
y = zeros(1, length(m));

for i = 1:length(m)
    for j = 1:length(n1)
        if m(i)==n1(j)
            newx1(i)=x1(j);
        end
    end
end

for i = 1:length(m)
    for j = 1:length(n2)
        if m(i)==n2(j)
            newx2(i)=x2(j);
        end
    end
end

y = newx1 + newx2;

%Signal Shifting
n2=[-3,-2,-1,0,1,2,3];
x2=[1,2,3,4,7,8,9];

m = n2+3;

subplot(2,1,1); stem(n2,x2); xlabel('n'); ylabel('x(n)'); axis([-5,5,1,12]), set(gca, 'xtick', [-5 -4 -3 -2 -1 0 1 2 3 4 5]);
subplot(2,1,2); stem(m,x2); xlabel('n'); ylabel('x(n-3)'); axis([-5,5,1,12]), set(gca, 'xtick', [-5 -4 -3 -2 -1 0 1 2 3 4 5]);

%===============================
%Convolution
x = [1,2,3,1];
h = [1,2,1,-1];

X = [x,zeros(1,length(h))];
H = [h,zeros(1,length(x))];

for n=1:length(x)+length(h)-1
    C(n) = 0;

    for k=1:length(x)
        if(n-k+1>0)
            C(n) = C(n) + X(k) * H(n-k+1);
        end
    end
    
end

nyb = nx(1)+nh(1); 
nye = nx(length(x)) + nh(length(h));
ny = [nybx:nye];
 

% Frequency analysis ploting
w = [0:1:500]*pi/500; % [0, pi] axis divided into 501 points.
X = exp(j*w)+2+3.*exp(-j*w)+4.*exp(-j*2*w)+5.*exp(-j*3*w);
magX = abs(X); angX = angle(X); realX = real(X); imagX = imag(X);
subplot(2,2,1); plot(w/pi,magX); grid
xlabel('frequency in pi units'); title('Magnitude Part'); ylabel('Magnitude')
subplot(2,2,3); plot(w/pi,angX); grid
xlabel('frequency in pi units'); title('Angle Part'); ylabel('Radians')
subplot(2,2,2); plot(w/pi,realX); grid
xlabel('frequency in pi units'); title('Real Part'); ylabel('Real')
subplot(2,2,4); plot(w/pi,imagX); grid
xlabel('frequency in pi units'); title('Imaginary Part'); ylabel('Imaginary')
        
        
        
    
