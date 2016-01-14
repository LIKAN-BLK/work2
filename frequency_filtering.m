function [X0, X1] = frequency_filtering(eegT, eegNT)
eegTf =  apply_bandpass_filter(eegT);
eegNTf = apply_bandpass_filter(eegNT);
for band =1:size(eegTf,4)
    [V, lam] = man_filt(eegTf(:,:,:,band), eegNTf(:,:,:,band),0.5);
    for i=1:size(eegTf,3)
        tmp = [];
        for j=1:1 tmp = cat(2,tmp,(eegTf(:,:,i,band) * V(:,j))'); end;
        X1(i,:,band) = tmp(1,:);
    end;
    for i=1:size(eegNTf,3)
        tmp = [];
        for j=1:1 tmp = cat(2,tmp, (eegNTf(:,:,i,band) * V(:,j))'); end;
        X0(i,:,band) = tmp(1,:);
    end;
end;
[X0,X1] = dimension_reduction(X0,X1);
end
function [ eegP] = apply_bandpass_filter(eeg)
%frequencyFiltering This function apples bandpass filter to eeg data.
%
%   Detailed explanation goes here

% f_min = 1;
% f_max = 4;                0.94, 0.47
% band_width = 25;
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

function [X0,X1] = dimension_reduction(X0,X1)
    Xtmp = cat(1,X0,X1);
    Xtotal=[];
    for i = 1:size(X0,3)
        Xtotal = cat(2,Xtotal,Xtmp(:,:,i));
    end;
    [coeff, score, latent] = pca(Xtotal);
    Xpca = Xtotal*coeff(:,1:10);
    X0 = Xpca(1:size(X0),:,:);
    X1 = Xpca(size(X0)+1:end,:,:);
end
