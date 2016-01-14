function [sens_tr, spec_tr, sens_tst, spec_tst] = parallel_classifying(eegT, eegNT,beta)
eegTf =  apply_bandpass_filter(eegT);
eegNTf = apply_bandpass_filter(eegNT);
N0 = size(eegNT, 3);
N1 = size(eegT, 3);
Y = [ones(1,N0) 2*ones(1,N1)]';
CV = cvpartition(Y, 'k', 5);
trIdx = CV.training(1);
tstIdx = CV.test(1);
Ytr = Y(trIdx, :);
Ytst= Y(tstIdx, :);
N0tr = sum(Ytr == 1);
N1tr = sum(Ytr == 2);
N0tst = sum(Ytst == 1);
N1tst = sum(Ytst == 2);
for band = 1:size(eegNTf,4)
	[V, lam] = man_filt(eegTf(:,:,:,band), eegNTf(:,:,:,band),beta);
	parfor i=1:size(eegTf,3)
		tmp = [];
		for j=1:1 tmp = cat(2,tmp,(eegTf(:,:,i,band) * V(:,j))'); end;
		X1(i,:) = tmp(1,:);
	end;
	parfor i=1:size(eegNTf,3)
		tmp = [];
		for j=1:1 tmp = cat(2,tmp, (eegNTf(:,:,i,band) * V(:,j))'); end;
		X0(i,:) = tmp(1,:);
	end;
	X = [X0;X1];
	
	Xtr = X(trIdx, :);
    Xtst = X(tstIdx, :);
	
	%train
	obj = train_shrinkage(Xtr, Ytr);
    W = obj.W;
	Q = X*W;
    Q0 = Q(Y == 1); %non target
    Q1 = Q(Y == 2);
	ths = Q + eps;
    ths = sort(ths);
	parfor k = 1:length(ths)                
        sens_all(k) = length(find(Q1 <= ths(k))) / N1;
        spec_all(k) = length(find(Q0 > ths(k))) / N0;    
    end;
	idx = find(spec_all >= 0.95, 1, 'last');
    th_opt= ths(idx);
	
    %train accuracy
    Q=[];
	Q=Xtr*W;
	tmp_Ypred_train = arrayfun(@(x) x < th_opt,Q);
    tmp_Ypred0_train = tmp_Ypred_train(Ytr == 1);
    tmp_Ypred1_train = tmp_Ypred_train(Ytr == 2);
	Ypred_train(band,:) = tmp_Ypred_train;
    band_sens_tr(band) = length(find(tmp_Ypred1_train == 1))/N1tr;
    band_spec_tr(band) = length(find(tmp_Ypred0_train == 0))/N0tr;

    
    %test
	Q=[];
	Q=Xtst*W;
	tmp_Ypred_test = arrayfun(@(x) x < th_opt,Q);
    tmp_Ypred0_test = tmp_Ypred_test(Ytst == 1);
    tmp_Ypred1_test = tmp_Ypred_test(Ytst == 2);
    Ypred_test(band,:) = tmp_Ypred_test;
    
    band_sens_tst(band) = length(find(tmp_Ypred1_test == 1))/N1tst;
    band_spec_tst(band) = length(find(tmp_Ypred0_test == 0))/N0tst;
end; 
test=1;


%Ypred_train=sum(tmp_Ypred_train,1);
% Ypred0_train = Ypred_train(Ytr == 1);
% Ypred1_train = Ypred_train(Ytr == 2);
% 
% sens_tr = length(find(Ypred1_train > size(eegNTf,4)/2))/N1tr;
% spec_tr  = length(find(Ypred0_train <= size(eegNTf,4)/2))/N0tr;
% 
% 
% Ypred_test=sum(tmp_Ypred_test,1);
% Ypred0_test = Ypred_test(Ytst == 1);
% Ypred1_test = Ypred_test(Ytst == 2);
% 
% sens_tst = length(find(Ypred1_test > size(eegNTf,4)/2))/N1tst;
% spec_tst = length(find(Ypred0_test <= size(eegNTf,4)/2))/N0tst;

end 

% function [sens_tr, spec_tr, sens_tst, spec_tst] = parallel_classifying(eegT, eegNT,beta)
% eegTf =  apply_bandpass_filter(eegT);
% eegNTf = apply_bandpass_filter(eegNT);
% N0 = size(eegNT, 3);
% N1 = size(eegT, 3);
% Y = [ones(1,N0) 2*ones(1,N1)]';
% CV = cvpartition(Y, 'k', 5);
% trIdx = CV.training(1);
% tstIdx = CV.test(1);
% Ytr = Y(trIdx, :);
% Ytst= Y(tstIdx, :);
% N0tr = sum(Ytr == 1);
% N1tr = sum(Ytr == 2);
% N0tst = sum(Ytst == 1);
% N1tst = sum(Ytst == 2);
% for band = 1:size(eegNTf,4)
% 	[V, lam] = man_filt(eegTf(:,:,:,band), eegNTf(:,:,:,band),beta);
% 	parfor i=1:size(eegTf,3)
% 		tmp = [];
% 		for j=1:1 tmp = cat(2,tmp,(eegTf(:,:,i,band) * V(:,j))'); end;
% 		X1(i,:) = tmp(1,:);
% 	end;
% 	parfor i=1:size(eegNTf,3)
% 		tmp = [];
% 		for j=1:1 tmp = cat(2,tmp, (eegNTf(:,:,i,band) * V(:,j))'); end;
% 		X0(i,:) = tmp(1,:);
% 	end;
% 	X = [X0;X1];
% 	
% 	Xtr = X(trIdx, :);
%     Xtst = X(tstIdx, :);
% 	
% 	%train
% 	obj = train_shrinkage(Xtr, Ytr);
%     W = obj.W;
% 	Q = X*W;
%     Q0 = Q(Y == 1); %non target
%     Q1 = Q(Y == 2);
% 	ths = Q + eps;
%     ths = sort(ths);
% 	parfor k = 1:length(ths)                
%         sens_all(k) = length(find(Q1 <= ths(k))) / N1;
%         spec_all(k) = length(find(Q0 > ths(k))) / N0;    
%     end;
% 	idx = find(spec_all >= 0.95, 1, 'last');
%     th_opt= ths(idx);
% 	
%     %train accuracy
%     Q=[];
% 	Q=Xtr*W;
% 	tmp_Ypred_train(band,:) = arrayfun(@(x) x < th_opt,Q);
%     
% 	%test
% 	Q=[];
% 	Q=Xtst*W;
% 	tmp_Ypred_test(band,:) = arrayfun(@(x) x < th_opt,Q);
% end; 
% 
% Ypred_train=sum(tmp_Ypred_train,1);
% Ypred0_train = Ypred_train(Ytr == 1);
% Ypred1_train = Ypred_train(Ytr == 2);
% 
% sens_tr = length(find(Ypred1_train > size(eegNTf,4)/2))/N1tr;
% spec_tr  = length(find(Ypred0_train <= size(eegNTf,4)/2))/N0tr;
% 
% 
% Ypred_test=sum(tmp_Ypred_test,1);
% Ypred0_test = Ypred_test(Ytst == 1);
% Ypred1_test = Ypred_test(Ytst == 2);
% 
% sens_tst = length(find(Ypred1_test > size(eegNTf,4)/2))/N1tst;
% spec_tst = length(find(Ypred0_test <= size(eegNTf,4)/2))/N0tst;
% 
% end

function [ eegP] = apply_bandpass_filter(eeg)
f_min = 1;
f_max = 25;               
band_width = 5;
f_order = 4;
fs=500;
for f_lo = f_min:f_max
    d = designfilt('bandpassiir', 'FilterOrder',2*f_order, 'HalfPowerFrequency1',f_lo, 'HalfPowerFrequency2',f_lo+band_width, 'SampleRate',fs);
    for ch = 1:size(eeg,2)
        parfor tr = 1:size(eeg,3)
            eegP(:,ch,tr,f_lo)=filter(d,eeg(:,ch,tr));
        end;
    end;
end; 
end
