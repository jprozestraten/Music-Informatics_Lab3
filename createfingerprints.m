%% Set constants
root = 'queries/';
file = 'Q1.wav';

tWindow = 0.050;
tHop = 0.010;

dt = 0.1;
nBands = 25;

%% Load query
[query,Fs] = audioread([root,file]);

%% Plot spectrogram
nWindow = floor(tWindow*Fs);
tOverlap = tWindow-tHop;
nOverlap = floor(tOverlap*Fs);

window = hann(nWindow,'periodic');

figure(1)
spectrogram(query,window,nOverlap,4*nWindow,Fs,'yaxis')
    ylim([0 5])

[S,f,t] = spectrogram(query,window,nOverlap,4*nWindow,Fs);

%% Find anchors
nIntervals = floor(t(end)/dt);

% Initialize matrices
fAnchors = zeros(nBands,nIntervals);
tAnchors = zeros(nBands,nIntervals);
valAnchors = zeros(nBands,nIntervals);

for i = 1:nBands
    for j = 1:nIntervals
        % Time and frequency range in which to look for max value
        fRange = (f >= (i-1)*Fs/2/nBands) & (f < i*Fs/2/nBands);
        tRange = (t >= (j-1)*dt) & (t < j*dt);
        
        % Bit of a work-around to get the proper indices from max
        fTmp = f(fRange);
        tTmp = t(tRange);
        
        % Get max value and indices
        [valMax,iiMax] = max(S(fRange,tRange),[],'all','linear');
        [iiF,iiT] = ind2sub([length(fTmp),length(tTmp)],iiMax);
        
        % Store anchor values
        fAnchors(i,j) = fTmp(iiF);
        tAnchors(i,j) = tTmp(iiT);
        valAnchors(i,j) = valMax;
    end
end
% Clear loop variables
clear i j fRange tRange fTmp tTmp valMax iiMax iiF iiT

% Convert values to decibel
valAnchors = 10*log10(abs(valAnchors));

%% Plot anchors in spectrogram
figure(1)
hold on
plot(tAnchors(:),fAnchors(:)*1e-3,'x','LineWidth',1)
hold off