% Quadratic Chirp Signal Parameters
a1 = 10;
a2 = 3;
a3 = 3;
A = 10;

% Instantaneous frequency after 1 sec is
maxFreq = a1 + 2 * a2 + 3 * a3;
% Nyquist frequency guess: 2 * max instantaneous frequency
nyqFreq = 2 * maxFreq;
% Sampling frequency
samplFreq = 5 * nyqFreq; 
samplIntrvl = 1 / samplFreq;

% Time samples
timeVec = 0:samplIntrvl:1.0;
% Number of samples
nSamples = length(timeVec);

% Generate the signal
sigVec = crcbgenqcsig(timeVec, A, [a1, a2, a3]);

% Plot the signal 
figure;
plot(timeVec, sigVec, 'Marker', '.', 'MarkerSize', 10, 'LineWidth', 1.5);
xlabel('Time (sec)');
ylabel('Amplitude');
title('Sampled Signal');
grid on;
saveas(gcf, 'sampled_signal.png');

% Plot the periodogram
% Length of data 
dataLen = timeVec(end) - timeVec(1);
% DFT sample corresponding to Nyquist frequency
kNyq = floor(nSamples / 2) + 1;
% Positive Fourier frequencies
posFreq = (0:(kNyq-1)) * (1 / dataLen);
% FFT of signal
fftSig = fft(sigVec);
% Discard negative frequencies
fftSig = fftSig(1:kNyq);

% Plot periodogram
figure;
plot(posFreq, abs(fftSig), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('|FFT|');
title('Periodogram');
grid on;
saveas(gcf, 'periodogram.png');

% Plot a spectrogram
winLen = 0.2; % sec
ovrlp = 0.1; % sec
% Convert to integer number of samples 
winLenSmpls = floor(winLen * samplFreq);
ovrlpSmpls = floor(ovrlp * samplFreq);
[S, F, T] = spectrogram(sigVec, winLenSmpls, ovrlpSmpls, [], samplFreq);

% Enhanced Spectrogram Plot
figure;
imagesc(T, F, abs(S)); axis xy;
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
title('Spectrogram');
colorbar;
set(gca, 'FontSize', 12);
saveas(gcf, 'spectrogram.png');

% Function to generate quadratic chirp signal
function sigVec = crcbgenqcsig(timeVec, A, params)
    % Generates a quadratic chirp signal
    % timeVec: vector of time stamps
    % A: amplitude of the signal
    % params: vector of chirp parameters [a1, a2, a3]
    
    a1 = params(1);
    a2 = params(2);
    a3 = params(3);
    
    % Compute the phase of the chirp
    phaseVec = a1 * timeVec + a2 * timeVec.^2 + a3 * timeVec.^3;
    
    % Generate the signal
    sigVec = A * cos(2 * pi * phaseVec);
end
