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

spectrogram(query,window,nOverlap,4*nWindow,Fs,'yaxis')
    ylim([0 5])
[S,f,t] = spectrogram(query,window,nOverlap,4*nWindow,Fs);

%% Find anchors
% for i = 
