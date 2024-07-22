sampFreq = 1024;
nSamples = 2048;
timeVec = (0:(nSamples - 1))/sampFreq;
a1 = 10;
a2 = 3;
a3 = 10;
A = 10;
siglen = (nSamples - 1)/sampFreq;
maxFreq = a1 + 2*a2*siglen + 3*a3*siglen^2;
disp(["the maximum frequency of the quadratic chirp is ", num2str(maxFreq)]);
sigVec = crcbgenqcsig(timeVec, A, [a1,a2,a3]);
filtorder = 30;
b = fir1(filtorder,(maxFreq/2)/(sampFreq/2));
filtsig = fftfilt(b, sigVec);
figure;
hold on;
plot(timeVec,sigVec);
plot(timeVec,filtsig);
xlabel('Time (sec)');
ylabel('Amplitude');
title('Sampled Signal');
function sigVec = crcbgenqcsig(dataX,snr,qcCoefs)
% Generate a quadratic chirp signal
% S = CRCBGENQSIG(X,SNR,C)
% Generates a quadratic chirp signal S. X is the vector of
% time stamps at which the samples of the signal are to be computed. SNR is
% the matched filtering signal-to-noise ratio of S and C is the vector of
% three coefficients [a1, a2, a3] that parametrize the phase of the signal:
% a1*t+a2*t^2+a3*t^3. 
%Soumya D. Mohanty, May 2018
phaseVec = qcCoefs(1)*dataX + qcCoefs(2)*dataX.^2 + qcCoefs(3)*dataX.^3;
sigVec = sin(2*pi*phaseVec);
sigVec = snr*sigVec/norm(sigVec);
end
