%%
close all;
clc;
clear;
%% ----------Reading Audio Signal (Hindi, contains 10 words; Counting 0 to 9)
[data, fs] = audioread('Voice sample\3.wav');
data = data./max(data); % Normalization
% plot(data,'b')
% title('Original Signal', 'fontsize',18)
%% ----------Removing silence part-------------------------------
data = silence_removal(data, fs)';

%% ------------Framing + Energy + MFCC----------------------------
%% ------------ Framing ------------------------------------------
f_duration = 0.025; % 25 msec
f_size = f_duration * fs;
n_frame = length(f_size/2+1: f_size/2: length(data))-2; % with 50% overlap
%% ------------- Energy ------------------------------------------
energy = zeros(1, n_frame); % Initialisation
%% ------------- preemphasis Filter--------------------------------
preemph = [1, -0.95];
%% ------------- Mel Filter Bank-----------------------------------
NFFT = 2^12;
K = NFFT/2+1;            % length of each filter
M = 26;                  % number of filters
hz2mel = @(hz)(1127*log(1+hz/700)); % Hertz to mel warping function
mel2hz = @(mel)(700*exp(mel/1127)-700); % mel to Hertz warping function
[ H1, freq ] = trifbank( M, K, [0 fs/2], fs, hz2mel, mel2hz );
% subplot(4,1,4)
% plot(freq, H1)
% title('26 Mel-spaced filter-bank', 'fontsize',18)
% axis tight
%% ------------- MFCC---------------------------------------------------
raw_c = zeros(n_frame, M); % For storing raw coefficients
j = 1;
for i = f_size/2+1: f_size/2: length(data)% With 50% Overlap
    if(i+f_size < length(data))
        frame = data(i-f_size/2:i+f_size/2); % taking datas inbetween the word "tin": 1000:1000+f_size
        energy(j) = sum(frame.^2); % Energy
        frame = filter(preemph, 1, frame); % Preemphasis Filter
        frame = frame.*hamming(length(frame)); % Hamming Window
        frame_f = abs(fft(data, NFFT));
        
        frame_f = frame_f(1:end/2+1);
        
        %  Periodogram estimate of the power spectrum
        frame_p = frame_f.^2./NFFT;
        
        % MFCC Filter + Log + Normalisation
        for k = 1:M
           raw_c(j,k) = log10(sum(frame_p'.*H1(k,:))/sum(H1(k,:))); % /sum(H1(k,:)) for normalising
        end
        
%         if(j == 250)
%             plot(raw_c(250,:),'b')
%             title('Filter-bank energy of 250th frame', 'fontsize',18)
%             axis tight
%         end
        
        % fprintf('%d \t %d\n',i-f_size/2, i+f_size/2)
        j = j + 1;
%         disp(j)
    end
end
% %%
%subplot(4,1,3)
plot(energy)
title('Energy of the signal', 'fontsize',18)
axis tight


%%
%-------------------DCT----------------------------
DCT_c = zeros(n_frame, M); % For storing DCT coefficients
for k = 1:n_frame
    j = 0;
    for i = 1:M
           DCT_c(k,i) = sum(raw_c(k,:).*cos(pi*j*((1:M)-0.5)./M))*sqrt(2/M);
           j = j+1;
    end
end
%% -------------------- MFCC Coefficients-----------------------------------
mfcc_coeff = zeros(n_frame, 12);
delta1_mfcc_coeff = zeros(n_frame-1, 12);
delta2_mfcc_coeff = zeros(n_frame-2, 12);

%--- Energy------------
delta1_energy = zeros(1,n_frame-1);
delta2_energy = zeros(1,n_frame-2);

mfcc_coeff(1,:) = DCT_c(1,1:12);
for i = 2:n_frame
    mfcc_coeff(i,:) = DCT_c(i,1:12);
    delta1_mfcc_coeff(i-1,:) = mfcc_coeff(i,:) - mfcc_coeff(i-1,:);
    delta1_energy(i-1) = energy(i) - energy(i-1);
end

plot(mfcc_coeff(250,:),'b')
title('Discrete Cosine Transform(DCT)', 'fontsize',18)
axis tight
for i=1:12
    txt = num2str(i);
    text(i,DCT_c(250,i)-0.5,txt, 'FontSize',10)
end

for i = 2:n_frame-1
    delta2_mfcc_coeff(i-1,:) = delta1_mfcc_coeff(i,:) - delta1_mfcc_coeff(i-1,:);
    delta2_energy(i-1) = delta1_energy(i) - delta1_energy(i-1);
end





