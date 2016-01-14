function [ufeats, usens] = eye_selectFeats(sens1, Nf, f_channels)
%
% select the best 1D-features
% ufeats: [t, channel]
%

% ------- select best features -------------
[best_sens, idx] = sort(sens1(:), 'descend');
nBest = min([2000 length(idx)]);
[feats(:,1), feats(:,2)] = ind2sub(size(sens1), idx(1:nBest));
feats = sortrows(feats, [2 1]);

nfeats = 0;
c = 1;
f = feats(1, :);
for i = 2:size(feats, 1)
    if (abs(feats(i, 1) - f(1)) < 5) && (feats(i, 2) == f(2))
        c = c + 1;
    elseif c >= 3
        nfeats = nfeats + 1;
        ufeats(nfeats, :) = feats(i - round(c/2), :);        
        usens(nfeats) = sens1(ufeats(nfeats, 1), ufeats(nfeats, 2));
        c = 1;
    end
    f = feats(i, :);
end

ufeats = ufeats(ufeats(:,1) > 60, :); % select only informative times
usens = usens(ufeats(:,1) > 60);
ufeats = ufeats(ufeats(:,2) <= 58, :); % select only informative channels
usens = usens(ufeats(:,2) <= 58);

% sift features
%f_channels = 1:58;
%f_channels = [8 11:20 22 23 24 39 40 42:53 55 56]; % P, O, C electrodes
idx = zeros(size(ufeats,1), 1);
for i = 1:length(f_channels)
    idx = idx | (ufeats(:,2) == f_channels(i));
end
ufeats = ufeats(idx, :);
usens = usens(idx);

[usens, idx] = sort(usens, 'descend');
ufeats = ufeats(idx, :);
Nf = min([size(ufeats, 1) Nf]);
ufeats = ufeats(1:Nf, :);
usens = usens(1:Nf);

%disp(['Channels: ' num2str(f_channels)])
