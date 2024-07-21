%% Initialize and Setup
% Add paths to folders containing signal and noise generation codes
addpath('../SIGNALS');
addpath('../NOISE');

%% Parameters for Data Realization
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1)) / sampFreq;

% Parameters for signal
a1 = 9.5;
a2 = 2.8;
a3 = 3.2;
A = 10;

%% Supply PSD Values
noisePSD = @(f) (f >= 100 & f <= 300) .* (f - 100) .* (300 - f) / 10000 + 1;
dataLen = nSamples / sampFreq;
kNyq = floor(nSamples / 2) + 1;
posFreq = (0:(kNyq-1)) * (1 / dataLen);
psdPosFreq = noisePSD(posFreq);

%% Generate Data Realization
% Generate noise and signal
noiseVec = statgaussnoisegen(nSamples, [posFreq(:), psdPosFreq(:)], 100, sampFreq);
sig4data = crcbgenqcsig(timeVec, A, [a1, a2, a3]);

% Normalize signal to SNR=10
[sig4data, ~] = normsig4psd(sig4data, sampFreq, psdPosFreq, 10);
dataVec = noiseVec + sig4data;

%% Plot Data and Spectrogram
figure;

% Time-domain plot
subplot(3,1,1);
plot(timeVec, dataVec);
hold on;
plot(timeVec, sig4data);
xlabel('Time (sec)');
ylabel('Data');
title('Data Realization for Calculation of GLRT');
legend('Data', 'Signal');

% Frequency-domain plot
subplot(3,1,2);
datFFT = abs(fft(dataVec));
sigFFT = abs(fft(sig4data));
plot(posFreq, datFFT(1:kNyq));
hold on;
plot(posFreq, sigFFT(1:kNyq));
xlabel('Frequency (Hz)');
ylabel('Periodogram');
legend('Data FFT', 'Signal FFT');

% Spectrogram
subplot(3,1,3);
[S, F, T] = spectrogram(dataVec, 64, 60, [], sampFreq);
imagesc(T, F, abs(S));
axis xy;
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
title('Spectrogram');

%% Compute GLRT
% Generate the template signal
sigVec = crcbgenqcsig(timeVec, A, [a1, a2, a3]);
[templateVec, ~] = normsig4psd(sigVec, sampFreq, psdPosFreq, 1);

% Calculate inner product of data with template
llr = innerprodpsd(dataVec, templateVec, sampFreq, psdPosFreq);

% GLRT is the square of the log-likelihood ratio
GLRT = llr^2;
disp(['GLRT Value: ', num2str(GLRT)]);



