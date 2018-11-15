%clc;
clear;
close();
%%
% Continuous periodic and random signal
A = 1;
w1 = 10;
w2 = 60;
t = 0:0.1/w2:2*pi;
theta = 0;

x1 = A*cos(w1*t + theta);
x2 = A*cos(w2*t + theta);%2*(rand(1,length(t))-0.5);

subplot(2,1,1)
plot(t,x1)
title('Low frequency: Continuous Signal \fontsize{46}')
xlabel('t = 0 - 2*Pi')
axis([0 2*pi -1 1])
subplot(2,1,2)
plot(t,x2)
title('High frequency: Continuous Signal \fontsize{46}')
xlabel('t = 0 - 2*Pi')
axis([0 2*pi -1 1])

%Continuous and discrete signal
A = 1;
f1 = 1;
fs = 2;

t = 0:0.01/(f1):2*pi;
n = 0:2*pi;
theta = 0;

x1 = A*cos(2*pi*f1*t + theta);
x2 = A*cos(2*pi*(f1/fs)*n + theta);

subplot(2,1,1)
plot(t,x1)
title('Continuous Signal')
xlabel('t = 0:2*Pi')
axis([0 2*pi -1 1])
subplot(2,1,2)
stem(n,x2)
title('Discrete Signal')
xlabel('t = 0:2*Pi')
axis([0 2*pi -1 1])

%% Impulse Signal --------------------------------------------------

A = 1;
f1 = 1;
%fs = 2;

t = 0:0.01/(f1):2*pi;
n = 0:2*pi;
theta = 0;

x1 = A*cos(2*pi*f1*t + theta);
x2 = A*cos(2*pi*n + theta);
x2 = ones(1,length(n)).*x2;

subplot(2,1,1)
plot(t,x1)
title('Continuous Signal')
xlabel('t = 0 - 2*Pi')
axis([0 2*pi -1 1])
subplot(2,1,2)
stem(n,x2)
title('Impulse Signal')
xlabel('t = 0 - 2*Pi')
axis([0 2*pi -1 1])

% 	
%% Calculate Power and to plot Power spectral density-------------------
syms t
x=input('Enter function in t: ','s');
T = input('Enter Limit :: T = ');   %Total evaluation time
Ts= 1; %0.001;   %Sampling time =&gt; 1000 samples per second
t1 = -T:0.01:T;
t=0:Ts:T; %define simulation time
x = eval(x);
power = sum(x)/(2*length(t1)+1);
energy = sum(x);
fprintf('Energy = %d\n', energy)
fprintf('Power = %d\n', power)

if(power<0.01)
  fprintf('It is Energy Signal.\n')
else
  fprintf('It is Power Signal.\n')
end



%% -------------------Convolution------------------------------------------
%---LInear Convolution--------------------------------------------

temp = magic(5);
a = temp(1,:);
b = temp(2,1:end);

convol = conv(a,b);
subplot(2,2,1)
stem(a)
title(['a = [',num2str(a),']'])

subplot(2,2,2)
stem(b)
title(['b = [',num2str(b),']'])

subplot(2,2,3)
stem(convol)
title(['Linear Convolution = [',num2str(convol),']'])

%% -----Circular Convolution: length of both the signal should be equal------------------------------------

a_pad = [a, zeros(1,abs(length(a)-length(b)))];
b_pad = [b, zeros(1,abs(length(a)-length(b)))];
c_convol = ifft(fft(a_pad).*fft(b_pad));
subplot(2,2,4)
stem(c_convol)
title(['Circular Convolution = [',num2str(c_convol),']'])

%% -----Circular Convolution equivalence to linear-----------------
a_cl = [a, zeros(1,length(a)+length(b)-1)];
b_cl = [b, zeros(1,length(a)+length(b)-1)];
cl_convol = ifft(fft(a_cl).*fft(b_cl));
subplot(2,3,5)
stem(cl_convol(1:(length(a)+length(b)-1)))
title(['cl_convol = [',num2str(cl_convol(1:(length(a)+length(b)-1))),']'])


%% ------------------------Correlation----------------------------
%------------Cross & Auto Correlation ----------------
a = rand(1,20);
b = a(1:5);
%b = rand(1,5);
correlation = xcorr(a(6:end),b);

subplot(3,4,1:4)
plot(a)
title('Signal A')

[max, i] = max(correlation);
subplot(3,4,ceil(i/10)+4)
plot(b)
title('Signal A[1:5]')

subplot(3,4,9:12)
stem(correlation)
title('Correlation Coefficient')

%% --------------------DFT fourier transform of speech signal | Spectrum-------------------

cd 'Voice Sample'
[data,fs] = audioread('segment3_1.wav');
cd ..
l = length(data);
NFFT = 1024;
k = 0:NFFT-1;
f = (fs/NFFT).*k;
xf = abs(fft(data, NFFT));

subplot(2,1,1);
plot(data);
title('Input Speech Signal');
subplot(2,1,2);
plot(f(1:end/2+1), xf(1:NFFT/2+1));
title('Single Sided Spectrum of the Speech Signal');
% 
%% ---------------------Spectrogram------------------
cd 'Voice Sample'
[data,fs] = audioread('segment3_1.wav');
cd ..
NFFT = 1024;
specgram(data, NFFT, fs)
% 
%% ----------Power spectral Analysis--------------

[data,fs] = audioread('Voice Sample/segment3_1.wav');
pwelch(data,[],[],[],fs)

%% ----- Addition of signals with different frequencies------------
clc;
close();
t = 0:0.01:5;
y1 = sin(2*pi*5*t);
subplot(3,1,1)
plot(t,y1)
title('High Frequency Signal')
xlabel('y1 = sin(2*pi*5*t)')

subplot(3,1,2)
y2 = sin(2*pi*t);
plot(t,y2)
title('Low Frequency Signal')
xlabel('y2 = sin(2*pi*t)')

subplot(3,1,3)
plot(t,y1+y2)
title('Addition of two signal')
xlabel('Y = Y1 + Y2')
