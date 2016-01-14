function res = eye_classify(eeg, params, fixDur, sRate, A1_ch, A2_ch)
%
% eeg: [Ntimes * Nchannels]
%

T = size(eeg, 1);
nChannels = size(eeg, 2);
W = params.W;
th = params.th;
ufeats = params.feats;

% preprocessing
eeg_p = eye_preprocess(eeg, fixDur, sRate, A1_ch, A2_ch);

X = zeros(1, size(ufeats, 1));
for i = 1:size(ufeats, 1);
    t = ufeats(i, 1);
    ch = ufeats(i, 2);
    X(i) = eeg_p(t, ch);    
end

Q = X*W;
res = Q < th; % 1: target, 0: non target

