% colGaussNoiseDemo.m
% Generate colored Gaussian noise and analyze its properties

% Parameters
Fs = 1000;            % Sampling frequency
N = 1024;             % Number of samples
factor_increase = [1, 2, 4]; % Factors to increase the number of samples

% Generate white Gaussian noise
wgn_noise = randn(N, 1);

% Filter design (example: low-pass filter)
fc = 100; % Cutoff frequency
[b, a] = butter(2, fc/(Fs/2)); % 2nd order Butterworth filter

for factor = factor_increase
    % Increase number of samples
    N_new = N * factor;
    wgn_noise = randn(N_new, 1);
    
    % Generate colored Gaussian noise
    colored_noise = filter(b, a, wgn_noise);
    
    % Estimate PSD using periodogram
    [Pxx, f] = periodogram(colored_noise, [], [], Fs);
    
    % Target PSD (for comparison, typically known or desired PSD)
    [H, w] = freqz(b, a, N_new, Fs);
    target_PSD = abs(H).^2;
    
    % Plot target and estimated PSD
    figure;
    subplot(2, 1, 1);
    plot(w, target_PSD, 'r', 'LineWidth', 2);
    hold on;
    plot(f, Pxx, 'b');
    title(['Target and Estimated PSD (N = ' num2str(N_new) ')']);
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
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
