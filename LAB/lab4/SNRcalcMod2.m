%% How to normalize a signal for a given SNR
% We will normalize a signal such that the Likelihood ratio (LR) test for it has
% a given signal-to-noise ratio (SNR) in noise with a given Power Spectral 
% Density (PSD). [We often shorten this statement to say: "Normalize the
% signal to have a given SNR." ]

%%
% Path to folder containing signal and noise generation codes
addpath ../SIGNALS
addpath ../NOISE

%%
% This is the target SNR for the LR
snr = 10;

%%
% Data generation parameters
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1))/sampFreq;

%%
% Generate your assigned signal (modify as needed)
% Example: a sine wave with frequency 5 Hz
freq = 5;
sigVec = sin(2*pi*freq*timeVec);

%%
% Read the initial LIGO design sensitivity PSD from file
psdData = load('../NOISE/iLIGOSensitivity.txt');
freqs = psdData(:,1);
psdVals = psdData(:,2);

%%
% Generate the PSD vector to be used in the normalization. Should be
% generated for all positive DFT frequencies.
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);

% Interpolate the PSD values to match the DFT frequencies
interpPsdVals = interp1(freqs, psdVals, posFreq, 'linear', 'extrap');

figure;
plot(posFreq, interpPsdVals);
axis([0,posFreq(end),0,max(interpPsdVals)]);
xlabel('Frequency (Hz)');
ylabel('PSD ((data unit)^2/Hz)');

%% Calculation of the norm
% Norm of signal squared is inner product of signal with itself
normSigSqrd = innerprodpsd(sigVec,sigVec,sampFreq,interpPsdVals);
% Normalize signal to specified SNR
sigVec = snr*sigVec/sqrt(normSigSqrd);

%% Test
% Obtain LLR values for multiple noise realizations
nH0Data = 1000;
llrH0 = zeros(1,nH0Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),interpPsdVals(:)],100,sampFreq);
    llrH0(lp) = innerprodpsd(noiseVec,sigVec,sampFreq,interpPsdVals);
end
% Obtain LLR for multiple data (=signal+noise) realizations
nH1Data = 1000;
llrH1 = zeros(1,nH1Data);
for lp = 1:nH0Data
    noiseVec = statgaussnoisegen(nSamples,[posFreq(:),interpPsdVals(:)],100,sampFreq);
    % Add normalized signal
    dataVec = noiseVec + sigVec;
    llrH1(lp) = innerprodpsd(dataVec,sigVec,sampFreq,interpPsdVals);
end
%%
% Signal to noise ratio estimate
estSNR = (mean(llrH1)-mean(llrH0))/std(llrH0);

figure;
histogram(llrH0);
hold on;
histogram(llrH1);
xlabel('LLR');
ylabel('Counts');
legend('H_0','H_1');
title(['Estimated SNR = ',num2str(estSNR)]);

%%
% A noise realization
figure;
plot(timeVec,noiseVec);
xlabel('Time (sec)');
ylabel('Noise');
%%
% A data realization
figure;
plot(timeVec,dataVec);
hold on;
plot(timeVec,sigVec);
xlabel('Time (sec)');
ylabel('Data');

%% Plot the periodogram of the signal
figure;
periodogram(sigVec,[],[],sampFreq);
title('Periodogram of the Signal');

%% Plot the periodogram of the data
figure;
periodogram(dataVec,[],[],sampFreq);
title('Periodogram of the Data');

%% Plot the spectrogram of the data
figure;
spectrogram(dataVec,128,120,128,sampFreq,'yaxis');
title('Spectrogram of the Data');



