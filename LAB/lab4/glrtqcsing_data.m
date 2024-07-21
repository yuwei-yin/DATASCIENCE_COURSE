% Main Script
clear; clc;

% Define parameters
a1 = 10;
a2 = 3;
a3 = 3;
sampFreq = 1024;
nSamples = 2048;
dataLen = nSamples / sampFreq;
posFreq = (0:(floor(nSamples / 2))) * (1 / dataLen);
noisePSD = @(f) (f >= 100 & f <= 300) .* (f - 100) .* (300 - f) / 10000 + 1;

% Load data realizations
dataFiles = {'DETEST/data1.txt', 'DETEST/data2.txt', 'DETEST/data3.txt'};
numRealizations = numel(dataFiles);
GLRTValues = zeros(numRealizations, 1);

for i = 1:numRealizations
    data = loadData(dataFiles{i});
    GLRTValues(i) = computeGLRT(data, (0:(nSamples-1))/sampFreq, 10, [a1, a2, a3], noisePSD(posFreq), sampFreq);
    
    % Plot data and signal
    figure;
    timeVec = (0:(nSamples-1))/sampFreq;
    sig4data = crcbgenqcsig(timeVec, 10, [a1, a2, a3]);
    [sig4data, ~] = normsig4psd(sig4data, sampFreq, noisePSD(posFreq), 10);
    
    plot(timeVec, data, 'b', 'DisplayName', 'Data');
    hold on;
    plot(timeVec, sig4data, 'r', 'DisplayName', 'Signal');
    xlabel('Time (sec)');
    ylabel('Amplitude');
    title(sprintf('Data and Signal for Realization %d', i));
    legend;
    grid on;
    
    % Plot frequency spectrum
    figure;
    dataFFT = abs(fft(data));
    sigFFT = abs(fft(sig4data));
    plot(posFreq, dataFFT(1:floor(nSamples/2)+1), 'b', 'DisplayName', 'Data FFT');
    hold on;
    plot(posFreq, sigFFT(1:floor(nSamples/2)+1), 'r', 'DisplayName', 'Signal FFT');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title(sprintf('Frequency Spectrum for Realization %d', i));
    legend;
    grid on;
    
    % Plot spectrogram
    figure;
    [S, F, T] = spectrogram(data, 64, 60, [], sampFreq);
    imagesc(T, F, abs(S)); axis xy;
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
    title(sprintf('Spectrogram for Realization %d', i));
    colorbar;
end

% Generate M realizations under H0
M = 1000; % Start with M = 1000 and increase as needed
GLRT_H0 = zeros(M, numRealizations);

for j = 1:M
    noiseVec = statgaussnoisegen(nSamples, [posFreq(:), noisePSD(posFreq(:))], 100, sampFreq);
    for i = 1:numRealizations
        GLRT_H0(j, i) = computeGLRT(noiseVec, (0:(nSamples-1))/sampFreq, 10, [a1, a2, a3], noisePSD(posFreq), sampFreq);
    end
end

% Estimate significance
significance = zeros(numRealizations, 1);
for i = 1:numRealizations
    significance(i) = mean(GLRT_H0(:, i) >= GLRTValues(i));
end

disp('GLRT Values:');
disp(GLRTValues);
disp('Significance Estimates:');
disp(significance);

% Profile the code (optional)
% profile on;
% Run the whole script or key sections here
% profile viewer;

%% Function Definitions
% Function to load data
function data = loadData(fileName)
    % Load data from a text file
    data = load(fileName);
end

% Function to compute GLRT
function GLRT = computeGLRT(dataVec, timeVec, A, params, psdPosFreq, sampFreq)
    % Generate the quadratic chirp signal
    sigVec = crcbgenqcsig(timeVec, A, params);
    
    % Normalize the signal to unit norm
    [templateVec, ~] = normsig4psd(sigVec, sampFreq, psdPosFreq, 1);
    
    % Calculate the inner product
    llr = innerprodpsd(dataVec, templateVec, sampFreq, psdPosFreq);
    
    % Ensure llr is a scalar or handle accordingly
    if ismatrix(llr)
        llr = llr(:); % Convert to a column vector if necessary
    end
    
    % Calculate GLRT (squared LLR)
    GLRT = sum(llr.^2); % Use sum to handle vector or matrix llr
end


