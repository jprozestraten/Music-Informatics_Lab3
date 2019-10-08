%% Compare fingerprint of query with itself
root = 'queries/';
file = 'Q1.wav';

t_0 = 2; % second;
dt = 1; % second;

% Read audiofile
[~,Fs] = audioread([root,file],[1 1]);
query = audioread([root,file],[floor(Fs*t_0) floor(Fs*(t_0+dt))]+1);

% Compute time-hash pairs of query
[tFingerprint,hFingerprint] = fingerprint(query,Fs);

% Compare fingerprint against itself
[iQuery,iDocument] = comparehash(hFingerprint,hFingerprint);

% Plot results
figure(1)
scatter(tFingerprint(iQuery),tFingerprint(iDocument),'kx')
    title('Comparison query with itself')
    xlabel('Time in query signal (s)')
    ylabel('Time in query signal (s)')
    xlim([0 max(tFingerprint)])
    ylim([0 max(tFingerprint)])

%% Compare fingerprint of query with original file
root = 'queries/';
file = 'Q1.wav';

t_0 = 2; % second;
dt = 1; % second;

% Read audiofile and compute time-hash pairs of query
[~,Fs] = audioread([root,file],[1 1]);
query = audioread([root,file],[floor(Fs*t_0) floor(Fs*(t_0+dt))]+1);

[tFingerprintQuery,hFingerprintQuery] = fingerprint(query,Fs);

% Read audiofile and compute time-hash pairs of document
[document,Fs] = audioread([root,file]);

[tFingerprintDocument,hFingerprintDocument] = fingerprint(document,Fs);

% Compare hashes against eachother
[iQuery,iDocument] = comparehash(hFingerprintQuery,hFingerprintDocument);

% Plot results
figure(2)
scatter(tFingerprintQuery(iQuery),tFingerprintDocument(iDocument),'kx')
    title('Comparison query with original file')
    xlabel('Time in query signal (s)')
    ylabel('Time in document signal (s)')
    xlim([0 max(tFingerprintQuery)])
    ylim([0 max(tFingerprintDocument)])

%% Compare fingeprint of query with another file
root = 'queries/';
file = 'Q1.wav';

t_0 = 2; % second;
dt = 1; % second;

% Read audiofile and compute time-hash pairs of query
[~,Fs] = audioread([root,file],[1 1]);
query = audioread([root,file],[floor(Fs*t_0) floor(Fs*(t_0+dt))]+1);

[tFingerprintQuery,hFingerprintQuery] = fingerprint(query,Fs);

% Read audiofile and compute time-hash pairs of different document
file = 'Q2.wav';
[document,Fs] = audioread([root,file]);

[tFingerprintDocument,hFingerprintDocument] = fingerprint(document,Fs);

% Compare hashes against eachother
[iQuery,iDocument] = comparehash(hFingerprintQuery,hFingerprintDocument);

% Plot results
figure(3)
scatter(tFingerprintQuery(iQuery),tFingerprintDocument(iDocument),'kx')
    title('Comparison query with different file')
    xlabel('Time in query signal (s)')
    ylabel('Time in document signal (s)')
    xlim([0 max(tFingerprintQuery)])
    ylim([0 max(tFingerprintDocument)])

%% 
function [iQuery,iDocument] = comparehash(query,document)
lenQuery = size(query,1);
lenDocument = size(document,1);

% Initialize empty vectors
iDocument = [];
iQuery = [];

% Loop through the query hashes
for i = 1:lenQuery
    % Subtract query hash from all document hashes, and sum the elements of
    % that result
    comp = sum(document-query(i,:),2); % Sum along 2nd dimension to sum hash elements
    comp = ~comp; % Returns all zeros
    
    iDocument = [iDocument;find(comp)]; % Append index vector with indices of nonzero elements
    iQuery = [iQuery;repmat(i,[sum(comp) 1])]; % Vector becomes [iQuery;i;i;...]. Amount of i's is the amount of matches
end
end