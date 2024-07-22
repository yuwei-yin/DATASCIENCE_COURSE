% crcbpso_script.m

% Step 1: Generate the data realization
nSamples = 512;
Fs = 512;
snr = 10;
a1 = 10;
a2 = 3;
a3 = 3;

% Create dataX
dataX = (0:(nSamples-1)) / Fs;

% Reset random number generator to generate the same noise realization,
% otherwise comment this line out
rng('default');

% Parameters for new function as a struct
P = struct('a1', a1, 'a2', a2, 'a3', a3);

% Generate data realization
[dataY, sig] = crcbgenqcsig_new(dataX, snr, P);

% Step 2: Run the modified crcbqcpso function
% Search range of phase coefficients
rmin = [1, 1, 1];
rmax = [15, 5, 5];

% Number of independent PSO runs
nRuns = 8;

% Input parameters for crcbqcpso
inParams = struct('dataX', dataX, ...
                  'dataY', dataY, ...
                  'dataXSq', dataX.^2, ...
                  'dataXCb', dataX.^3, ...
                  'rmin', rmin, ...
                  'rmax', rmax);

% Run PSO
outStruct = crcbqcpso(inParams, struct('maxSteps', 2000), nRuns);

% Step 3: Make the figure
figure;
hold on;

% Plot the data realization
plot(dataX, dataY, '.', 'DisplayName', 'Data');

% Plot the true signal
plot(dataX, sig, 'DisplayName', 'True Signal');

% Plot the best estimated signal from each PSO run
colors = lines(nRuns);
for lpruns = 1:nRuns
    plot(dataX, outStruct.allRunsOutput(lpruns).estSig, 'Color', colors(lpruns, :), 'LineWidth', 1.5, 'DisplayName', ['Run ', num2str(lpruns)]);
end

% Plot the best estimated signal
plot(dataX, outStruct.bestSig, 'Color', 'k', 'LineWidth', 2.0, 'DisplayName', 'Best Estimated Signal');

% Add legend
legend('show');

% Add titles and labels
title('PSO Signal Estimation');
xlabel('Time (s)');
ylabel('Amplitude');

% Display estimated parameters
disp(['Estimated parameters: a1=', num2str(outStruct.bestQcCoefs(1)), ...
                             '; a2=', num2str(outStruct.bestQcCoefs(2)), ...
                             '; a3=', num2str(outStruct.bestQcCoefs(3))]);

hold off;

