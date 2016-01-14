function [spec1, sens1] = eye_1Dfeats(eegT, eegNT)
%
% calc classification accuray in 1D-feature spaces
%

nChannels = size(eegT, 2) - 1;
T = size(eegT, 1);
N0 = size(eegNT, 3);
N1 = size(eegT, 3);

% ------------- calc classification accuracy ---------------------
spec1 = zeros(T, nChannels);
sens1 = zeros(T, nChannels);
for t = 1:T
    for ch = 1:nChannels
        X0 = squeeze(eegNT(t, ch, :));
        X1 = squeeze(eegT(t, ch, :));
        
        X = [X0; X1];
        Y = [ones(1,N0) 2*ones(1,N1)]';
        obj = train_shrinkage(X,Y);
        W = obj.W;

        Q0 = W'*X0;
        Q1 = W'*X1;
        ths = [Q0; Q1] + eps;
        ths = sort(ths);
        for k = 2:length(ths)                
            sens(k) = length(find(Q1 <= ths(k))) / N1;
            spec(k) = length(find(Q0 > ths(k))) / N0;
            if (spec(k) < 0.95)
                idx = k - 1;
                break;
            end
        end;
        %idx = find(spec >= 0.95, 1, 'last');
        spec1(t, ch) = spec(idx);
        sens1(t, ch) = sens(idx);
        %[acc(t, ch), idx] = max((sens * N1 + spec * N0) / (N1 + N0));
        th = ths(idx);
    end
end    

spec1;

