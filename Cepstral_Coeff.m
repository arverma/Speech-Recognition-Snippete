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
subplot(2,2,1)
frame_f = abs(frame_f);
plot(xfft, frame_f(1:nfft/2))
legend('Fourier Transform');
%%
%----------Ceptral Transform Coefficient
subplot(2,2,2)
plot(xfft, frame_f(1:nfft/2))
hold on;
frame_f = log10(frame_f); %---LOG
plot(xfft,abs(frame_f(1:nfft/2)))
legend('Fourier Transform', 'Log');
%%
%---------IDFT
frame_f = ifft(abs(frame_f(1:nfft/2)),nfft);
subplot(2,2,3)
plot(abs(frame_f(1:nfft/2)/max(frame_f)))
legend('IDFT');
%%
%---------IDFT
subplot(2,2,4)
stem(abs(frame_f(1:13)))
legend('Ceptral Coefficient');


