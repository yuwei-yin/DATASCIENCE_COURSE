% Load the data
data = load('testData.txt');
times = data(:, 1);
values = data(:, 2);

% Plot the time series data
figure;
plot(times, values);
hold on;
xline(5.0, 'r--', 'Signal starts');
title('Time Series Data');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original Data');
hold off;
% Use data before 5.0 seconds
noise_values = values(times < 5.0);

% Estimate the noise PSD
[psd, frequencies] = pwelch(noise_values);

% Plot the PSD
figure;
semilogy(frequencies, psd);
title('Noise Power Spectral Density');
xlabel('Frequency (Hz)');
ylabel('PSD (V^2/Hz)');
grid on;
% Design the whitening filter
whitening_filter = 1 ./ sqrt(psd);

% Plot the whitening filter response
figure;
plot(frequencies, 20 * log10(abs(whitening_filter)));
title('Whitening Filter Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
grid on;
% Apply the whitening filter
whitened_values = filter(whitening_filter, 1, values);

% Plot the time series before and after whitening
figure;
subplot(2, 1, 1);
plot(times, values);
hold on;
xline(5.0, 'r--', 'Signal starts');
title('Time Series Data Before Whitening');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original Data');
hold off;

subplot(2, 1, 2);
plot(times, whitened_values, 'g');
hold on;
xline(5.0, 'r--', 'Signal starts');
title('Time Series Data After Whitening');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Whitened Data');
hold off;
% Compute the spectrograms
[~, f_original, t_original, Sxx_original] = spectrogram(values, 256, 250, 256, 1/(times(2) - times(1)));
[~, f_whitened, t_whitened, Sxx_whitened] = spectrogram(whitened_values, 256, 250, 256, 1/(times(2) - times(1)));

% Plot the spectrograms
figure;
subplot(2, 1, 1);
imagesc(t_original, f_original, 10*log10(abs(Sxx_original)));
axis xy;
title('Spectrogram Before Whitening');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

subplot(2, 1, 2);
imagesc(t_whitened, f_whitened, 10*log10(abs(Sxx_whitened)));
axis xy;
title('Spectrogram After Whitening');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
