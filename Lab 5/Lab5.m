
%% Sampling and reconstruction demo
clear,clc,close all;

%% Parameters
F = 30;     % frequency of signal [Hz]
Fs = 2*F;   % sampling rate [Hz]
Ts = 1/Fs;  % sampling period [sec]

%% Generate "continuous time" signal and discrete time signal
tc = 0:1e-4:5/F;        % CT axis
xc = cos(2*pi*F*tc);    % CT signal
td = 0:Ts:5/F;          % DT axis
xd = cos(2*pi*F*td);    % DT signal
N = length(td);         % number of samples

%% Reconstruction by using the formula:
% xr(t) = sum over n=0,...,N-1: x(nT)*sin(pi*(t-nT)/T)/(pi*(t-nT)/T)
% Note that sin(pi*(t-nT)/T)/(pi*(t-nT)/T) = sinc((t-nT)/T)
% sinc(x) = sin(pi*x)/(pi*x) according to MATLAB
xr = zeros(size(tc));
sinc_train = zeros(N,length(tc));
for t = 1:length(tc)
    for n = 0:N-1
        sinc_train(n+1,:) = sin(pi*(tc-n*Ts)/Ts)./(pi*(tc-n*Ts)/Ts);
        xr(t) = xr(t) + xd(n+1)*sin(pi*(tc(t)-n*Ts)/Ts)/(pi*(tc(t)-n*Ts)/Ts);
        
    end
end

%% Plot the results
figure
hold on
grid on
plot(tc,xc)
stem(td,xd)
plot(tc,xr, 'r')
xlabel('Time [sec]')
ylabel('Amplitude')

%% Sinc train visualization    
figure
hold on
grid on
plot(tc,xd'.*sinc_train)
stem(td,xd)
plot(tc,xc,'r', 'LineWidth',2)
xlabel('Time [sec]')
ylabel('Amplitude')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%program to find the DFT/IDFT of a sequence without using the inbuilt functions
close all;
clear all;
N=6;
n=0:N-1;
xn=cos((pi*n)/3);
ck=zeros(1,N);

for k=0:N-1
    for n=0:N-1
        ck(k+1)=ck(k+1)+(xn(n+1)*exp((-i)*2*pi*k*n/N));
    end
end

ck=ck./N;

% Signal plot
Nn=0:N-1;
subplot(221);
stem(Nn,xn);
ylabel ('Amplitude');
xlabel ('Time Index');

% Power Spectrum
subplot(222);
stem(Nn,ck.*ck);
ylabel ('ck');
xlabel ('N');

% Find the magnitudes of individual DFT points
k=0:1:N-1;
magnitude=abs(ck);
subplot(223);
stem(k,magnitude);
ylabel ('Amplitude');
xlabel ('K');

% Find the phases of individual DFT points
phase=angle(ck);
subplot(224);
stem(k,phase);
ylabel ('Phase');
xlabel ('K');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Homework
%1. Discrete Fourier Series (Synthesis equation)
%2. Discrete Fourier Transform (DFT)
%3. Inverse Discrete Fourier Transform (IDFT)
%4. DFT as linear combination
%5. Circular Convolution

