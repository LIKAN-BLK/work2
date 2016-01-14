function [ output_args ] = statistics( eegT, eegNT,beta, ufeats, A1_ch, A2_ch, fixationDuration, sRate,chans_for_clf )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i = 1:size(eegT,1)-2
    eegTp = eegT(i:i+2) - repmat(mean(eegT(i:i+2),1),3,1)


end

