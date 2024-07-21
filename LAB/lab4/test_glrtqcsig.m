% Define the parameters
aVec = [10, 3, 3];
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1))/sampFreq;

% Generate some example data
sigVec = crcbgenqcsig(timeVec, 1, aVec);
noisePSD = @(f) (f>=100 & f<=300).*(f-100).*(300-f)/10000 + 1;
dataLen = nSamples/sampFreq;
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdPosFreq = noisePSD(posFreq);
noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdPosFreq(:)],100,sampFreq);
dataVec = noiseVec + sigVec;

% Calculate the GLRT value
glrtVal = glrtqcsig(dataVec, timeVec, psdPosFreq, aVec);
disp(['GLRT value: ', num2str(glrtVal)]);
', num2str(glrtVal)]);
test_glrtqcsig
