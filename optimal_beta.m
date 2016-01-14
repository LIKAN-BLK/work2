function [sensivity,specificity] = optimal_beta( eegT, eegNT, ufeats, A1_ch, A2_ch, fixationDuration, sRate, chans_for_clf )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
max_sens_test(1) = 0;
max_sens_tr(1) = 0;
max_spec_test(1) = 0;
max_spec_tr(1) = 0;
tmpmax_sens_test(1) = 0;
tmpmax_sens_tr(1) = 0;
tmpmax_spec_test(1) = 0;
tmpmax_spec_tr(1) = 0;
for dt = 0:(size(eegT,1)-1)
    for tau=1:(size(eegT,1)-dt)    
        for i=1:11
%             [params, spec, sens, acc] = eye_train1(eegT(tau:tau+dt,:,:), eegNT(tau:tau+dt,:,:),0.5, ufeats, A1_ch, A2_ch, fixationDuration, sRate,chans_for_clf);
%             sens_tst(tau,dt+1,:) = sens.tst;
%             sens_tr(tau,dt+1,:) = sens.tr;
%             spec_tst(tau,dt+1,:) = spec.tst;
%             spec_tr(tau,dt+1,:) = spec.tr;
            
            [params, spec, sens, acc] = eye_train1(eegT(tau:tau+dt,:,:), eegNT(tau:tau+dt,:,:),(i-1)/10, ufeats, A1_ch, A2_ch, fixationDuration, sRate,chans_for_clf);
            sens_tst3D(tau,dt+1,i)=sens.tst(1);
            sens_tr3D(tau,dt+1,i)=sens.tr(1);
            spec_tst3D(tau,dt+1,i)=spec.tst(1);
            spec_tr3D(tau,dt+1,i)=spec.tr(1);
            if(tmpmax_sens_test(1) <= sens.tst(1))
                tmpmax_sens_test = sens.tst;
                tmpsens_test_beta = (i-1)/10;
            end;
            if(tmpmax_sens_tr(1) <= sens.tr(1))
                tmpmax_sens_tr = sens.tr;
                tmpsens_tr_beta = (i-1)/10;
            end;
            if(tmpmax_spec_test(1) <= spec.tst(1))
                tmpmax_spec_test = spec.tst;
                tmpspec_test_beta = (i-1)/10;
            end;
            if(tmpmax_spec_tr(1) <= spec.tr(1))
                tmpmax_spec_tr = spec.tr;
                tmpspec_tr_beta = (i-1)/10;
            end;
        end;
        
        if(max_sens_test(1) < tmpmax_sens_test(1))
            max_sens_test = tmpmax_sens_test;
            max_sens_test_tau = tau; max_sens_test_dt = dt; sens_test_beta = tmpsens_test_beta;
        end;
        if(max_sens_tr(1) < tmpmax_sens_tr(1))
            max_sens_tr = tmpmax_sens_tr;
            max_sens_tr_tau = tau; max_sens_tr_dt = dt; sens_tr_beta = tmpsens_tr_beta;
        end;
        if(max_spec_test(1) < tmpmax_spec_test(1))
            max_spec_test = tmpmax_spec_test;
            max_spec_test_tau = tau; max_spec_test_dt = dt; spec_test_beta = tmpspec_test_beta;
        end;
        if(max_spec_tr(1) < tmpmax_spec_tr(1))
            max_spec_tr = tmpmax_spec_tr;
            max_spec_tr_tau = tau; max_spec_tr_dt = dt; spec_tr_beta = tmpspec_tr_beta;
        end;
    end;
end;
for dt = 0:(size(eegT,1)-1)
    for tau=1:(size(eegT,1)-dt)
        [tmp_max, tmp_index] = max(sens_tst3D(tau,dt+1,:));
        image_sens_tst(tau,dt+1) = (tmp_index-1)/10;
        [tmp_max, tmp_index] = max(sens_tr3D(tau,dt+1,:));
        image_sens_tr(tau,dt+1) = (tmp_index-1)/10;
        [tmp_max, tmp_index] = max(spec_tst3D(tau,dt+1,:));
        image_spec_tst(tau,dt+1) = (tmp_index-1)/10;
        [tmp_max, tmp_index] = max(spec_tr3D(tau,dt+1,:));
        image_spec_tr(tau,dt+1) = (tmp_index-1)/10;
    end;
end;

%imagesc(sens_tst), xlabel('dt'), ylabel('Tau'), title('beta'),grid minor, set(gca, 'xtick', 1:1:size(eegT,1)),set(gca, 'ytick', 1:1:size(eegT,1)) ;
sensivity.test.value = max_sens_test;
sensivity.test.tau = max_sens_test_tau;
sensivity.test.dt = max_sens_test_dt;
sensivity.test.beta = sens_test_beta;

sensivity.train.value = max_sens_tr;
sensivity.train.tau = max_sens_tr_tau;
sensivity.train.dt = max_sens_tr_dt;
sensivity.train.beta = sens_tr_beta;

specificity.test.value = max_spec_test;
specificity.test.tau = max_spec_test_tau;
specificity.test.dt = max_spec_test_dt;
specificity.test.beta = spec_test_beta;

specificity.train.value = max_spec_tr;
specificity.train.tau = max_spec_tr_tau;
specificity.train.dt = max_spec_tr_dt;
specificity.train.beta = spec_tr_beta;

% figure(1),imagesc(image_sens_tst), xlabel('dt'), ylabel('Tau'), title('sens tst beta'),grid minor, set(gca, 'xtick', 1:1:size(eegT,1)),set(gca, 'ytick', 1:1:size(eegT,1)) ;
% figure(2),imagesc(image_sens_tr), xlabel('dt'), ylabel('Tau'), title('sens tr beta'),grid minor, set(gca, 'xtick', 1:1:size(eegT,1)),set(gca, 'ytick', 1:1:size(eegT,1)) ;
% figure(3),imagesc(image_spec_tst), xlabel('dt'), ylabel('Tau'), title('spec tst beta'),grid minor, set(gca, 'xtick', 1:1:size(eegT,1)),set(gca, 'ytick', 1:1:size(eegT,1)) ;figure(4),imagesc(image_spec_tr), xlabel('dt'), ylabel('Tau'), title('spec tr beta'),grid minor, set(gca, 'xtick', 1:1:size(eegT,1)),set(gca, 'ytick', 1:1:size(eegT,1)) ;
% figure(4),imagesc(image_spec_tr), xlabel('dt'), ylabel('Tau'), title('spec tr beta'),grid minor, set(gca, 'xtick', 1:1:size(eegT,1)),set(gca, 'ytick', 1:1:size(eegT,1)) ;figure(4),imagesc(image_spec_tr), xlabel('dt'), ylabel('Tau'), title('spec tr beta'),grid minor, set(gca, 'xtick', 1:1:size(eegT,1)),set(gca, 'ytick', 1:1:size(eegT,1)) ;

end

