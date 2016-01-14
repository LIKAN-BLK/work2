function [eegTwp,eegNTwp] = bandpass_filter(eegTp, eegNTp, sRate)
    lowFq = 3;
    hightFq = 30;
    num_bands=4;
    step = (hightFq-lowFq + 1)/num_bands; 
    for ch = 1:size(eegTp,2)
        for band_number = 1:num_bands
            for tr = 1:size(eegTp,3)
                eegTwp(:,(ch-1)*num_bands + band_number,tr) = windowed_bandpass_filter(eegTp(:,ch,tr), lowFq+(band_number-1)*step, lowFq+band_number*step-1, sRate);
            end;
            for tr = 1:size(eegNTp,3)
                eegNTwp(:,(ch-1)*num_bands + band_number,tr) = windowed_bandpass_filter(eegNTp(:,ch,tr), lowFq+(band_number-1)*step, lowFq+band_number*step-1, sRate);
            end;
        end;
    end;
%     Tzeros_cols=[];
%     Tzeros_rows=[];
%     for tr=1:size(eegTwp,3)
%         Tzeros_cols = union(Tzeros_cols, find(~sum(eegTwp(:,:,tr),1)));
%         Tzeros_rows = union(Tzeros_rows, find(~sum(eegTwp(:,:,tr),2)));
%     end;
%   
%     NTzeros_cols=[];
%     NTzeros_rows=[];
%     for tr=1:size(eegNTwp,3)
%         NTzeros_cols = union(NTzeros_cols,find(~sum(eegNTwp(:,:,tr),1)));
%         NTzeros_rows = union(NTzeros_rows,find(~sum(eegNTwp(:,:,tr),2)));
%     end;
%     
%     zeros_cols = union(Tzeros_cols,NTzeros_cols);
%     zeros_rows = union(Tzeros_rows,NTzeros_rows);
%     eegTwp(:,Tzeros_cols,:) = [];
%     eegTwp(Tzeros_rows,:,:) = [];
%     
%     eegNTwp(:, NTzeros_cols,:) = [];
%     eegNTwp(NTzeros_rows, :, :) = [];
%     
end


function eegWP = windowed_bandpass_filter(eeg, Fpass1, Fpass2, sRate)
    halfRate = sRate/2;
    [b,a] = butter(5,[Fpass1,Fpass2]/halfRate, 'bandpass'); 
    eegWP = filter(b,a,eeg);
    
end
