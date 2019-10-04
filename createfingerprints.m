%% Set constants
root = 'queries/';
file = 'Q1.wav';

% STFT settings
tWindow = 0.050; % seconds
tHop = 0.010;    % seconds

% Time-frequency box dimensions of anchors
dt = 0.1;    % seconds
nBands = 25; % amount of bands to divide spectrum into

% Relative targetzone
tTargetzone = [+0.100 +0.500]; % Addition-wise
fTargetzone = [2^-0.5 2^0.5];  % Multiplication-wise

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

%% Create fingerprint
fingerprint = cell(numel(tAnchors),2);

fingerprint(:,1) = num2cell(tAnchors(:));

hash = zeros(1,3);

% Loops through all anchors
for i = 1:numel(tAnchors)
    % Time and frequency range of anchors within targetzone
    fRange = (fAnchors >= fAnchors(i)*fTargetzone(1)) & (fAnchors < fAnchors(i)*fTargetzone(2));
    tRange = (tAnchors >= tAnchors(i)+tTargetzone(1)) & (tAnchors < tAnchors(i)+tTargetzone(2));
    
    % Range of anchors where both time and frequency are in targetzone
    range = fRange & tRange;
    
    % 
    
end