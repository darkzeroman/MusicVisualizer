function out =  musicvisualizer(file, useGPU, seconds)
if nargin < 2, useGPU = 1; end
if nargin < 3, seconds = 0; end
if seconds == 0
    [y, Fs] = wavread(file);      % data= y, sampling freq = Fs
else
    [y, Fs] = wavread(file,[1 44100*seconds]);   % for a certain range
end

% Common Variables
player = audioplayer(y,Fs);
y = y(:,2)*10;
if (useGPU)
    y = gsingle(y);
end
lenY = length(y);
lengthSong = lenY/Fs;
timeToJump = 2.2;
timeSkip = 350/44100;
goodTill = 0;

% Limits for each subplot
g2max=1;
g3max=13;
g4max=2;

% For the spheres
[xs,ys,zs] = sphere(10);
sphTrans = [0,0,0;10,-10,10;10,-10,-10;0,0,-15;-10,10,-10;-10,10,10];
xs = xs.*ones(size(xs),class(y));
ys = ys.*ones(size(ys),class(y));
zs = zs.*ones(size(zs),class(y));

if (useGPU)
    % GPU Graphing works a bit different, so slightly different points are
    % needed
    sphTrans = gsingle([0,0,0; -10,-10,0;  -10,0,0; 0,-15,0;0,0,-10; 0,-10,-10]);
    
end

% Variables to use in the loop
lastTimeFFT = 0;
time = 0;

val = [];
play(player);
myt = tic;
while time<lengthSong
    tic
    time = toc(myt);
    
    % For the FFT Analysis, which happens only once every 2.2 seconds
    if (or((goodTill <= time),(goodTill-time) < .2)) % Testing if we're past the last analyzed point, or close to it
        if (lengthSong - time > timeToJump)
            % There are more than 2.2 seconds left in the song
            startFrameFFT = ceil(time*44100);
            endFrameFFT = ceil((time+timeToJump)*44100);
            sigBandsFFT = analyze(y(startFrameFFT:endFrameFFT));
            goodTill = time + timeToJump;
            lastTimeFFT = time;
        else
            % If less then 2.2 seconds left, grab the last 2.2 seconds
            sigBandsFFT = analyze(y((ceil((lengthSong - timeToJump)*44100):ceil((lengthSong)*44100))));
            goodTill = lenY;
            lastTimeFFT = time;
        end
        
        % For the second plot, normalized individual bands
        upperlimit = ones(1,1,class(y))*g2max;
        if (useGPU)
            gsubplot(2,2,3);
            % upperlimit = gones(1, 1)*g2max;
            plot(upperlimit); % Plotting the upper limit
            ghold_on;
        else
            subplot(2,2,2);
            %upperlimit = ones(1, 1)*g2max;
            plot(upperlimit);
            hold on;
        end
        
        plot(upperlimit*0); % Plotting the lower limit
        i = [];
        for i=1:6
            c = abs(sigBandsFFT(:,i));
            c = c/max(c);
            plot(c);
        end
        
        if (useGPU)
            ghold_off
        else
            hold off;
        end
    end
    
    % waveformRange is for the wave form ranges, which is updated on every loop run
    % If at the beginning of song, just grabbing the first interval
    if (time < timeSkip)
        time = timeSkip+.001;
    end
    waveformRange = [ceil((time-timeSkip)*44100) ceil((time+timeSkip)*44100)];
    
    % Making sure to not have a limit beyond song length
    if (waveformRange(2) >= lenY)
        waveformRange(2)= lenY;
    end
    waveformRange = waveformRange(1):waveformRange(end);
    
    % If the input range is empty, we're past the end of the song
    if isempty(waveformRange)
        break;
    end
    
    % First plot, waveform
    if (useGPU)
        gsubplot(2,2,1);
    else
        subplot(2,2,1);
    end
    plot(y(waveformRange)/max(abs(y(waveformRange)))) %plotting the data with limits -1 to 1
    
    time = toc(myt);
    
    %%% For the third plot, six spheres
    normSigBandsFFT = zeros(size(sigBandsFFT),class(y));
    
    if (useGPU)
        gsubplot(2,2,2);
    else
        subplot(2,2,3);
        
    end
    % Find the difference between last FFT and current time to find the
    % frame to use
    deltatime = time-lastTimeFFT;
    deltaFrame = ceil(deltatime*44100);
    if deltaFrame == 0 %frame index has to be atleast 1
        deltaFrame = 1;
    end
    
    sigBandsFFT = abs(sigBandsFFT);
    maxValuesFFT = max(sigBandsFFT);
    minValuesFFT = min(sigBandsFFT);
    
    %Normalizing with respect to the max values and initializing the
    %scatter plot variables
    X = ones(0, class(y));
    Y = ones(0,class(y));
    Z = ones(0,class(y));
    if (useGPU)
        gfor i=1:6
        %Signal Bands FFT, normalized (somewhat)
        normSigBandsFFT(:,i) =  2*sigBandsFFT(:,i)/maxValuesFFT(i);
        gend
        
    else
        for i=1:6
            normSigBandsFFT(:,i) =  2*sigBandsFFT(:,i)/maxValuesFFT(i);
        end
        
    end
    
    for i = 1:6
        multFactor = normSigBandsFFT(deltaFrame,i);
        X = [X; xs(:)*multFactor-sphTrans(i,1)];
        Y = [Y; ys(:)*multFactor-sphTrans(i,2)];
        Z = [Z; zs(:)*multFactor-sphTrans(i,3)];
    end
    
    X = [X;-1*g3max;1*g3max]; %limits
    Y = [Y;-1*g3max;1*g3max];
    Z = [Z;-1*g3max;1*g3max];
    scatter3(X(:),Y(:),Z(:))
    
    
    currNormSigBandsFFT = zeros(1,6,class(y));
    % Fourth plot, a single way to see the beat
    if (useGPU)
        gsubplot(2,2,4)
    else
        subplot(2,2,4)
    end
    
    for i=1:6
        currNormSigBandsFFT(i) =currNormSigBandsFFT(i)+ normSigBandsFFT(deltaFrame,i) ;
    end
    
    [a b] = max(maxValuesFFT);
    [c d] = min(minValuesFFT);
    
    
    currNormSigBandsFFT = currNormSigBandsFFT - .5*min(currNormSigBandsFFT);
    currNormSigBandsFFT = currNormSigBandsFFT/max(currNormSigBandsFFT);
    
    X = [xs(:)*currNormSigBandsFFT(b)]; %middle
    Y = [ys(:)*currNormSigBandsFFT(b)];
    Z = [zs(:)*currNormSigBandsFFT(b)];
    
    X = [X;-1*g4max;1*g4max]; %limits
    Y = [Y;-1*g4max;1*g4max];
    Z = [Z;-1*g4max;1*g4max];
    scatter3(X(:),Y(:),Z(:))
    
    
    pause(.01);
    %toc;
    val = [val toc];
    mean(val) % uncomment to calculating the average time in the loop.
    
end
out = 0;