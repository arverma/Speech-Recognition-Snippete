%% CEPSTRAL TRANSFORM COEFFICIENT
close all;
clc;
clear;
%% ----------Reading Audio Signal (Hindi, contains 10 words; Counting 0 to 9)
[data, fs] = audioread('Voice sample\3.wav');
data = data(1:floor(length(data)/25)); % word "tin": 3
data = data./max(data); % Normalization
subplot(2,3,1:3)
plot(data)
%%
%-----------Framing the audio signal
f_duration = 0.025;
f_size = f_duration *fs;
frame = data(1000:1000+f_size); % taking datas inbetween the word "tin": 1000:1000+f_size
hold on;
plot(1000:1000+f_size, data(1000:1000+f_size),'+r'); % To point out the frame in original signal
legend('Original', 'Frame')

subplot(2,3,4)
plot(frame);
legend('Original Frame')
%%
%-----------Preemphasis Filter, all zero filter
preemph = [1, -0.95];
frame = filter(preemph, 1, frame);
subplot(2,3,5)
plot(frame);
legend('Pre-emphasis Frame')
%%
%-----------Hamming Windowing
subplot(2,3,6)
frame = frame.*hamming(length(frame));
plot(frame,'g');
hold on;
plot(hamming(length(frame)));
legend('After Windowing', 'Hamming Window');
%%
%-----------DFT
nfft = 2.^nextpow2(length(frame));
xfft = fs.*(0:nfft/2-1)/nfft;
frame_f = fft(frame, nfft);
% frame_f2 = frame_f(1:nfft/2); 
figure;
subplot(3,2,1)
frame_f = abs(frame_f);
plot(xfft, frame_f(1:nfft/2))
legend('Fourier Transform');
%%
%----------Ceptral Transform Coefficient
subplot(3,2,2)
plot(xfft, frame_f(1:nfft/2))
hold on;
frame_f = log10(frame_f); %---LOG
plot(xfft,abs(frame_f(1:nfft/2)))
legend('Fourier Transform', 'Log');
%%
%---------Low Pass Filter
subplot(3,2,3)
order = 32;
fc = 1e3/fs/2;
h = fir1(order, fc);%---Low Pass Filter Impulse Response
hf = fft(h,nfft);
% hf = hf(1:nfft/2);
plot(xfft, abs(hf(1:nfft/2)))
legend('Low Pass filter');
mul = hf'.*frame_f;
subplot(3,2,4)
plot(xfft, abs(hf(1:nfft/2)))
hold on;
plot(xfft, abs(mul(1:nfft/2)))
legend('Low Pass filter', 'Filtered Signal');

%%
%---------IDFT
frame_f = ifft(mul,nfft);
subplot(3,2,5)
plot(abs(frame_f(1:nfft/2)/max(frame_f)))
legend('IDFT');
%%
%---------IDFT
subplot(3,2,6)
stem(abs(frame_f(1:13)))
legend('Ceptral Coefficient');


