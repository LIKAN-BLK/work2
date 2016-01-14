function optimal_time_window(eegT, eegNT, ufeats, A1_ch, A2_ch, fixationDuration, sRate, chans_for_clf  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
max_sens_test = 0;
tmpmax_sens_test = 0;
for dt = 0:(size(eegT,1)-1)
    for tau=1:(size(eegT,1)-dt)    
        
            [params, spec, sens, acc] = eye_train1(eegT(tau:tau+dt,:,:), eegNT(tau:tau+dt,:,:),0.5, ufeats, A1_ch, A2_ch, fixationDuration, sRate,chans_for_clf);
            sens_tst(tau,dt+1) = sens.tst(1);
            sens_tr(tau,dt+1) = sens.tr(1);
            spec_tst(tau,dt+1) = spec.tst(1);
            spec_tr(tau,dt+1) = spec.tr(1);
       
    end;
end;
[tau_sens_tst, dt_sens_tst ] =find(sens_tst == max(sens_tst(:)));
[tau_sens_tr,dt_sens_tr]= find(sens_tr == max(sens_tr(:)));
[tau_spec_tst, dt_spec_tst] = find(spec_tst == max(spec_tst(:)));
[tau_spec_tr,dt_spec_tr] = find(spec_tr == max(spec_tr(:)));
imagesc(sens_tst), xlabel('dt'), ylabel('Tau'), title('sens tst'),grid minor, set(gca, 'xtick', 1:1:size(eegT,1)),set(gca, 'ytick', 1:1:size(eegT,1)) ;
figure(2),imagesc(sens_tr), xlabel('dt'), ylabel('Tau'), title('sens tr'),grid minor, set(gca, 'xtick', 1:1:size(eegT,1)),set(gca, 'ytick', 1:1:size(eegT,1)) ;
figure(3),imagesc(spec_tst), xlabel('dt'), ylabel('Tau'), title('spec tst'),grid minor, set(gca, 'xtick', 1:1:size(eegT,1)),set(gca, 'ytick', 1:1:size(eegT,1)) ;
figure(4),imagesc(spec_tr), xlabel('dt'), ylabel('Tau'), title('spec tr'),grid minor, set(gca, 'xtick', 1:1:size(eegT,1)),set(gca, 'ytick', 1:1:size(eegT,1)) ;

end
