function eeg = eye_preprocess(eeg, fixDur, sRate, A1_ch, A2_ch)
%
% eeg: [Ntimes * Nchannels]
%

T = size(eeg, 1);
nChannels = size(eeg, 2);

% montage
ref = (eeg(:, A2_ch) + eeg(:, A1_ch)) / 2;
e1 = eeg - repmat(ref, 1, nChannels);    
e1(:, A2_ch) = -ref; %as Cz
eeg = e1;

% centering
% base = mean(eeg, 1);
% s = std(eeg, [], 1);
% %e1 = eeg - repmat(base, T, 1);
% e1 = (eeg - repmat(base, T, 1)) ./ repmat(s, T, 1);
% eeg = e1;

%smoothing
dec_n = 5;
samples_in_epoch = fixDur/1000*sRate;
e1 = zeros(samples_in_epoch/dec_n, nChannels);
for i = 1:nChannels
   e1(:, i) = decimate(eeg(:, i), dec_n);
end
eeg = e1;


% filtering
% [b, a] = butter(4, [0.1*2/fs 4*2/fs]);
% d = designfilt('bandpassfir', 'FilterOrder',20, 'CutoffFrequency1',0.1, 'CutoffFrequency2',5, 'SampleRate',fs);
