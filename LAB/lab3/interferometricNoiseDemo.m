% interferometricNoiseDemo.m
% Simulate noise for an interferometric detector using the initial LIGO sensitivity curve

% Parameters
Fs = 4096;            % Sampling frequency
N = 2^18;             % Number of samples (factor can be adjusted)
factor_increase = [1, 2, 4]; % Factors to increase the number of samples

% Obtain initial LIGO sensitivity curve data
% Frequency (Hz) and corresponding PSD values (1/Hz)
% Data can be obtained from LIGO documentation or online sources
f_data = [10, 20, 30, 40, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000]; % Example frequencies
psd_data = [1e-21, 1e-22, 5e-23, 2e-23, 1e-23, 5e-24, 2e-24, 1.5e-24, 1.2e-24, 1e-24, 9e-25, 8e-25, 7e-25, 6e-25, 5e-25, 4e-25, 3e-25, 2e-25]; % Example PSD

% Interpolate the PSD data for the desired frequency resolution
f_interpolated = linspace(min(f_data), Fs/2, N/2+1);
psd_interpolated = interp1(f_data, psd_data, f_interpolated, 'pchip');

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
