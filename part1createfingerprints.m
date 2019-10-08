%% Set constants
root = 'queries/';
file = 'Q1.wav';

% STFT settings
tWindow = 0.050; % seconds
tHop = 0.010;    % seconds

% Time-frequency box dimensions of anchors
dt = 0.100;    % seconds
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
    title(['STFT of ',file])
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
    title(['STFT of ',file,' with anchors'])
hold off

%% Create fingerprint
% Declares a cell array for the times of all anchors
% Will later be transformed to a list of hashes with times
allHashes = cell(numel(tAnchors),1);

% Number of hashes per anchor
nHashes = zeros(numel(tAnchors),1);

% Loops through all anchors
for i = 1:numel(tAnchors)
    % Time and frequency range of anchors within targetzone
    fRange = (fAnchors >= fAnchors(i)*fTargetzone(1)) & (fAnchors < fAnchors(i)*fTargetzone(2));
    tRange = (tAnchors >= tAnchors(i)+tTargetzone(1)) & (tAnchors < tAnchors(i)+tTargetzone(2));
    
    % Range of anchors where both time and frequency are in targetzone
    range = fRange & tRange;
    nHashes(i) = sum(range,'all');
    
    % Time and frequency of all valid anchors
    fTmp = fAnchors(range);
    tTmp = tAnchors(range);
    
    % Make hashes
    hashes = zeros(nHashes(i),3);
    
    hashes(:,1) = fAnchors(i);
    hashes(:,2) = fTmp;
    hashes(:,3) = tTmp - tAnchors(i);
    
    % Store hashes in cell
    allHashes{i} = hashes;
end

% Store data in vectors/matrices
tFingerprint = repelem(tAnchors(:),nHashes);
hFingerprint = cell2mat(allHashes);

% Clear loop variables
clear allHashes hashes fTmp tTmp fRange tRange range
