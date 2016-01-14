function [params, spec, sens, acc] = eye_train1(eegT, eegNT, beta, ufeats, fixationDuration, sRate, chans_for_clf)
%
%
%


[V, lam] = man_filt(eegT, eegNT,beta);


% for i=1:size(eegT,3)
%     tmp = [];
%     for j=1:1 tmp = cat(2,tmp,(eegT(:,:,i) * V(:,j))'); end;
%     X1(i,:) = tmp(1,:);
% end;
% 
% for i=1:size(eegNT,3)
%     tmp = [];
%     for j=1:1 tmp = cat(2,tmp, (eegNT(:,:,i) * V(:,j))'); end;
%     X0(i,:) = tmp(1,:);
% end;
[X0, X1] = frequency_filtering(eegT, eegNT);

nChannels = size(eegT, 2) - 1;
T = size(eegT, 1);
N0 = size(eegNT, 3);
N1 = size(eegT, 3);

%for i = 1:size(ufeats, 1);
%    %%freq = ufeats(i, 1);
%    t = ufeats(i, 1);
%    ch = ufeats(i, 2);
%    X0(:, i) = eegNT(t, ch, :);
%    X1(:, i) = eegT(t, ch, :);
%end

X = [X0; X1];
Y = [ones(1,N0) 2*ones(1,N1)]';
CV = cvpartition(Y, 'k', 5);
Wplot=[];
for i = 1:CV.NumTestSets
    trIdx = CV.training(i);
    tstIdx = CV.test(i);
    Xtr = X(trIdx, :);
    Xtst = X(tstIdx, :);
    Ytr = Y(trIdx, :);
    Ytst = Y(tstIdx, :);
    N0tr = sum(Ytr == 1);
    N1tr = sum(Ytr == 2);
    N0tst = sum(Ytst == 1);
    N1tst = sum(Ytst == 2);

    % train
    obj = train_shrinkage(Xtr, Ytr);
    W(:,:,i) = obj.W;
    Wplot(:,i) = obj.W;
    % calc threshold using all sample
    Q = X*W(:,:,i);
    Q0 = Q(Y == 1); %non target
    Q1 = Q(Y == 2); %target
    ths = Q + eps;
    ths = sort(ths);
    for k = 1:length(ths)                
        sens_all(k) = length(find(Q1 <= ths(k))) / N1;
        spec_all(k) = length(find(Q0 > ths(k))) / N0;    
    end;
    idx = find(spec_all >= 0.95, 1, 'last');
    th_opt(i) = ths(idx);
    
    % calc acc on train sample
    Q = Xtr*W(:,:,i);    
    Q0 = Q(Ytr == 1);
    Q1 = Q(Ytr == 2);
    sens_tr(i) = length(find(Q1 <= th_opt(i))) / N1tr;
    spec_tr(i) = length(find(Q0 > th_opt(i))) / N0tr;    
    acc_tr(i) = (sens_tr(i) * N1tr + spec_tr(i) * N0tr) / (N1tr + N0tr);    
    
%     ths = Q + eps;
%     ths = sort(ths);
%     for k = 1:length(ths)                
%         sens_tr(k) = length(find(Q1 <= ths(k))) / N1tr;
%         spec_tr(k) = length(find(Q0 > ths(k))) / N0tr;    
%     end;
%     %[acc_tr(freq, t, ch), idx] = max((sens_tr * N1tr + spec_tr * N0tr) / (N1tr + N0tr));
%     idx = find(spec_tr >= 0.97, 1, 'last');
%     spec_tr_opt(i) = spec_tr(idx);
%     sens_tr_opt(i) = sens_tr(idx);
%     acc_tr_opt(i) = (sens_tr_opt(i) * N1tst + spec_tr_opt(i) * N0tst) / (N1tst + N0tst);
%     th_opt = ths(idx);
    
    % test
    Q = Xtst*W(:,:,i);    
    Q0 = Q(Ytst == 1);
    Q1 = Q(Ytst == 2);
    sens_tst(i) = length(find(Q1 <= th_opt(i))) / N1tst;
    spec_tst(i) = length(find(Q0 > th_opt(i))) / N0tst;    
    acc_tst(i) = (sens_tst(i) * N1tst + spec_tst(i) * N0tst) / (N1tst + N0tst);
    %plot(1-spec_tr, sens_tr);    
end
nchannels = size(chans_for_clf,2);
chans_for_clf{1,strcmp(chans_for_clf,'A2')} = 'Cz';

spec.tr = [mean(spec_tr) std(spec_tr)];
sens.tr = [mean(sens_tr) std(sens_tr)];
acc.tr = [mean(acc_tr) std(acc_tr)];
spec.tst = [mean(spec_tst) std(spec_tst)];
sens.tst = [mean(sens_tst) std(sens_tst)];
acc.tst = [mean(acc_tst) std(acc_tst)];

params.W = mean(W, 3);
params.th = mean(th_opt);
params.feats = ufeats;
params.fixationDuration = fixationDuration;
params.sRate = sRate;

% fid = fopen('../res/acc3.txt', 'w');
% fprintf(fid, 'Train:\n');
% fprintf(fid, ' Sensitivity: %f +- %f\n', sens.tr);
% fprintf(fid, ' Specificity: %f +- %f\n', spec.tr);
% fprintf(fid, '\n');
% fprintf(fid, 'Test:\n');
% fprintf(fid, ' Sensitivity: %f +- %f\n', sens.tst);
% fprintf(fid, ' Specificity: %f +- %f\n', spec.tst);
% fclose(fid);