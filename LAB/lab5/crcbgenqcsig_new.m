function [dataY, sig] = crcbgenqcsig_new(dataX, snr, P)
    % Generate quadratic chirp signal using structure parameters
    sig = P.a1*dataX + P.a2*dataX.^2 + P.a3*dataX.^3;
    % Add noise
    noise = randn(size(dataX));
    dataY = sig + noise * std(sig) / 10^(snr / 20);
end
