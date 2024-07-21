% interferometricNoiseDemo.m
% Simulate noise for an interferometric detector using the initial LIGO sensitivity curve

% Load sensitivity curve data from the file
data = load('iLIGOSensitivity.txt');
frequencies = data(:, 1);
psd_values = data(:, 2);

% Parameters
Fs = 4096;              % Sampling frequency
N = 2^18;               % Initial number of samples
factor_increase = [1, 2, 4]; % Factors to increase the number of samples

% Interpolate the PSD data for the desired frequency resolution
f_interpolated = linspace(min(frequencies), Fs/2, N/2+1);
psd_interpolated = interp1(frequencies, psd_values, f_interpolated, 'pchip');

for factor = factor_increase
    % Increase number of samples
    N_new = N * factor;
    
    % Generate white Gaussian noise
    wgn_noise = randn(N_new, 1);
    
    % FFT of white noise
    WGN_fft = fft(wgn_noise);
    
    % Apply the interpolated PSD
    target_PSD = psd_interpolated .^ 0.5; % Square root to apply to amplitude spectrum
    target_PSD = [target_PSD, fliplr(target_PSD(2:end-1))]; % Mirror for negative frequencies
    colored_noise_fft = WGN_fft .* target_PSD(:);
    
    % Inverse FFT to get time-domain signal
    colored_noise = ifft(colored_noise_fft, 'symmetric');
    
    % Estimate PSD using periodogram
    [Pxx, f] = periodogram(colored_noise, [], [], Fs);
    
    % Plot target and estimated PSD
    figure;
    subplot(2, 1, 1);
    loglog(f_interpolated, psd_interpolated, 'r', 'LineWidth', 2);
    hold on;
    loglog(f, Pxx, 'b');
    title(['Target and Estimated PSD (N = ' num2str(N_new) ')']);
    xlabel('Frequency (Hz)');
    ylabel('Power Spectral Density (1/Hz)');
    legend('Target PSD', 'Estimated PSD');
    hold off;
    
    % Plot noise time series
    subplot(2, 1, 2);
    plot(colored_noise);
    title(['Colored Gaussian Noise Time Series (N = ' num2str(N_new) ')']);
    xlabel('Sample');
    ylabel('Amplitude');
    xlim([1, min(1000, N_new)]); % Zoom in for better visualization
    
    % Plot histogram
    figure;
    histogram(colored_noise, 'Normalization', 'pdf');
    hold on;
    x = linspace(min(colored_noise), max(colored_noise), 100);
    plot(x, normpdf(x, 0, std(colored_noise)), 'r', 'LineWidth', 2);
    title(['Histogram of Colored Gaussian Noise (N = ' num2str(N_new) ')']);
    xlabel('Amplitude');
    ylabel('Probability Density');
    legend('Histogram', 'Normal PDF');
    hold off;
end

% Comments on observations
disp('The estimated PSD becomes smoother with an increase in the number of samples due to the reduction in the variance of the periodogram estimate.');
disp('Examine the noise time series by zooming in. It should look more structured compared to WGN due to the applied filter.');
disp('The histogram should still resemble a normal distribution if the noise remains Gaussian after filtering.');

