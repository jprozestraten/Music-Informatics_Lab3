root = 'queries/';
file = 'Q1.wav';

t_0 = 2; % second;
dt = 1; % second;

[~,Fs] = audioread([root,file],[1 1]);
query = audioread([root,file],[floor(Fs*t_0) floor(Fs*(t_0+dt))]+1);

document = audioread([root,file]);

[tFingerprint,hFingerprint] = createfingerprint(query,Fs);

[iQuery,iDocument] = comparehash(hFingerprint,hFingerprint);

scatter(tFingerprint(iQuery),tFingerprint(iDocument),'kx')

function [iQuery,iDocument] = comparehash(query,document)
    lenQuery = size(query,1);
    lenDocument = size(document,1);
    
    iDocument = [];
    iQuery = [];
    for i = 1:lenQuery
        comp = sum(document-query(i,:),2);
        comp = ~comp;
        
        iDocument = [iDocument;find(comp)];
        iQuery = [iQuery;repmat(i,[sum(comp) 1])];
    end
end