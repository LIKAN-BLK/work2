function [eegTwp,eegNTwp] = windowed_preprocess(eegT, eegNT, window_size,f_channels)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    window_size = window_size-1;
    %[eegTwp,eegNTwp] = windowed_preprocess_mean(eegT, eegNT, window_size ); %(wsize:3, dec:5) - 0.55, %(wsize:5, dec:5) - 0.57
    
    %[eegTwp,eegNTwp] = windowed_preprocess_unbiased_varianse(eegT, eegNT, window_size ); %(wsize:3, dec:5) - 0.42; %(wsize:5, dec:5) -0.42
    %[eegTwp,eegNTwp] = windowed_preprocess_biased_varianse(eegT, eegNT, window_size ); %(wsize:3, dec:5) -0.42
    %[eegTwp,eegNTwp] = windowed_preprocess_correlation(eegT, eegNT, window_size,f_channels ); % (wsize:3, dec:5) - 0.62 (!!!!!)
    [eegTwp,eegNTwp] = windowed_preprocess_normalisation(eegT, eegNT, window_size ); %(wsize:3, dec:5) - 0.35,%(wsize:3, dec:5) -0.35
end


function [eegTwp,eegNTwp] = windowed_preprocess_mean(eegT, eegNT, window_size )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    for i = 1:(size(eegT,1)-window_size)
        eegTwp(i,:,:)=mean(eegT(i:i+window_size,:,:),1);
        eegNTwp(i,:,:)=mean(eegNT(i:i+window_size,:,:),1);
    end;

end


function [eegTwp,eegNTwp] = windowed_preprocess_unbiased_varianse(eegT, eegNT, window_size )
    for i = 1:(size(eegT,1)-window_size)
        eegTwp(i,:,:)=var(eegT(i:i+window_size,:,:),0,1);     %Не работает
        eegNTwp(i,:,:)=var(eegNT(i:i+window_size,:,:),0,1);
    end;
end
function [eegTwp,eegNTwp] = windowed_preprocess_biased_varianse(eegT, eegNT, window_size )
    for i = 1:(size(eegT,1)-window_size)
        eegTwp(i,:,:)=var(eegT(i:i+window_size,:,:),1,1);
        eegNTwp(i,:,:)=var(eegNT(i:i+window_size,:,:),1,1);
    end;
end


function [eegTwp,eegNTwp] = windowed_preprocess_correlation(eegT, eegNT, window_size,f_channels )
    f_channels = sort(f_channels);
    eegT = eegT(:,f_channels,:);
    eegNT = eegNT(:,f_channels,:);
    mask = triu(ones(size(eegT,2), size(eegT,2)),1);
        for it = 1:(size(eegT,1)-window_size)
            for tr = 1:size(eegT,3)
                tmp = corrcoef(eegT(it:it+window_size,:,tr));
                eegTwp(it,:,tr) = tmp(mask == 1);
            end;
            for tr = 1:size(eegNT,3)
                tmp = corrcoef(eegNT(it:it+window_size,:,tr));
                eegNTwp(it,:,tr) = tmp(mask == 1);
            end;    
        end;
    
    %end;
end
function [eegTwp,eegNTwp] = windowed_preprocess_normalisation(eegT, eegNT, window_size )
   for i = 1:(size(eegT,1)-window_size)
        eegTwp(i,:,:)=(eegT(i+floor(window_size/2),:,:)-mean(eegT(i:i+window_size,:,:),1))./var(eegT(i:i+window_size,:,:),1,1);
        eegNTwp(i,:,:)=(eegNT(i+floor(window_size/2),:,:)-mean(eegNT(i:i+window_size,:,:),1))./var(eegNT(i:i+window_size,:,:),1,1);
   end; 
end

