function output = analyze(y)
%y = gsingle(y(:,1)); %When using Timeit, uncomment this to force GPU mode

% Common Variables
maxSigFreq=44100;
hannWinLength = .4;
hannlen = hannWinLength*2*maxSigFreq;

% Finding out if we're using GPU or CPU depending on the input
if (isa(y,'garray'))
    usingGPU = true;   
else
    usingGPU = false;
end

bandRanges = zeros(12,1,class(y));
bandlimits=[0 200 400 800 1600 3200].*ones(1,6,class(y));
numBands = length(bandlimits);
yFFTed = fft(y);
lenY = length(yFFTed);
output = zeros(lenY,numBands,class(y)); %for storing the outputs


% Using band limits to find the ranges
if (usingGPU)
    i = [];
    gfor i = 1:numBands-1
        bandRanges(2*(i)-1) = floor(bandlimits(i)/maxSigFreq*lenY/2)+1;
        bandRanges(2*i) = floor(bandlimits(i+1)/maxSigFreq*lenY/2);
    gend
else
    for i = 1:numBands-1
        bandRanges(2*(i)-1) = floor(bandlimits(i)/maxSigFreq*lenY/2)+1;
        bandRanges(2*i) = floor(bandlimits(i+1)/maxSigFreq*lenY/2);
    end
end
% The end cases of the band limits are special cases
bandRanges(2*numBands-1) = floor(bandlimits(numBands)/maxSigFreq*lenY/2)+1;
bandRanges(numBands*2) = floor(lenY/2);


% Using the frequency bands to format the output
for i = 1:numBands
    start = bandRanges(2*i-1);
    stop = bandRanges(2*i);
    output(start:stop,i) = yFFTed(start:stop);
    output(lenY+1-stop:lenY+1-start,i) = yFFTed(lenY+1-stop:lenY+1-start);
end
output(1,1)=0;


sigTime = ones(size(output),class(y));    % signal in time domain
sigFreq = ones(size(sigTime),class(y));    % signal in frequency domain


% For the Hann
hann = zeros(lenY,1);
% For the (relatively) small value of hannlen, it's slower to use the GPU.
for a = 1:hannlen
    hann(a) = (cos(a*pi/hannlen/2)).^2;
end


% Going from frequency domain and then back with absolute values in the time domain
if (usingGPU)
    % Now converting to GSingle if needed
    hann = gsingle(hann);
    k = [];
    gfor k = 1:numBands
    sigFreq(:,k) = fft(abs(real(ifft(output(:,k)))));
    gend
else
    for k = 1:numBands
        sigFreq(:,k) = fft(abs(real(ifft(output(:,k)))));
    end
end

% Half-Hanning FFT * Signal FFT. And then back to time domain
if (usingGPU)
    i = [];
    gfor i = 1:numBands
    output(:,i) = real(ifft(sigFreq(:,i).*fft(hann)));
    gend
else
    for i = 1:numBands
        output(:,i) = real(ifft(sigFreq(:,i).*fft(hann)));
    end
end
