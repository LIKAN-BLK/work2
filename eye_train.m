function params = eye_train(eegT, eegNT, chans_for_clf, parameters, sRate, path, epoch_size)
%
% eegT, eegNT: [Ntimes, Nchannels, Ntrials]
%
A1_ch = read_param_struct({'A1'}, parameters);
A2_ch = read_param_struct({'A2'}, parameters);

% preprocess
for i = 1:size(eegT, 3)
    eegTp(:, :, i) = eye_preprocess(eegT(:, :, i), epoch_size, sRate, A1_ch, A2_ch);
end

for i = 1:size(eegNT, 3)
    eegNTp(:, :, i) = eye_preprocess(eegNT(:, :, i), epoch_size, sRate, A1_ch, A2_ch);
end

% select features
% -- calc 1D accuracies
[spec1, sens1] = eye_1Dfeats(eegTp, eegNTp);
% %figure; imagesc(sens1'); colorbar;
Nf = 110;
f_channels =  read_param_struct(chans_for_clf, parameters);

%[ufeats, usens] = eye_selectFeats(sens1, Nf, f_channels);


% -- user defined features
f_times = 1:size(eegTp,1);
%f_channels = 1:58;
%f_channels = [8 11:20 22 23 24 39 40 42:53 55 56]; % P, O, C electrodes
[tmp1, tmp2] = meshgrid(f_times, f_channels);
ufeats = [tmp1(:) tmp2(:)];

%frequency_filtering(eegTp, eegNTp)

% [eegTwp,eegNTwp] = windowed_preprocess(eegTp, eegNTp, 3,f_channels);
% [senswp,senswp] = optimal_beta( eegTwp, eegNTwp, ufeats, A1_ch, A2_ch, fixationDuration, sRate, chans_for_clf );

%[ tauwp,dtwp,betawp,max_sens_testwp ] = optimal_beta( eegTwp, eegNTwp, ufeats, A1_ch, A2_ch, fixationDuration, sRate, chans_for_clf );

%imagesc(sens_tst), xlabel('dt'), ylabel('Tau'), title('sens test'),grid minor, set(gca, 'xtick', 1:1:size(eegTp,1)),set(gca, 'ytick', 1:1:size(eegTp,1)) ;
%figure (2), imagesc(sens_tr), xlabel('dt'), ylabel('Tau'), title('sens train'), grid minor, set(gca, 'xtick', 1:1:size(eegTp,1)),set(gca, 'ytick', 1:1:size(eegTp,1)) ;
%figure (3), imagesc(spec_tst), xlabel('dt'), ylabel('Tau'), title('spec test'), grid minor, set(gca, 'xtick', 1:1:size(eegTp,1)),set(gca, 'ytick', 1:1:size(eegTp,1)) ;
%figure (4), imagesc(spec_tr), xlabel('dt'), ylabel('Tau'), title('spec train'), grid minor, set(gca, 'xtick', 1:1:size(eegTp,1)),set(gca, 'ytick', 1:1:size(eegTp,1)) ;

[params, spec, sens, acc] = eye_train1(eegTp, eegNTp,0.5, ufeats, sRate,chans_for_clf);
%[sens_tr, spec_tr, sens_tst, spec_tst] = parallel_classifying_perceptron(eegTp,eegNTp,0.5);
% disp('Test:');
% disp(sprintf(' Specificity: ', spec_tr));
% disp(sprintf(' Sensitivity: \n', sens_tr));
% disp('Test:');
% disp(sprintf(' Specificity: %f +- %f', spec_tst));
% disp(sprintf(' Sensitivity: %f +- %f\n', sens_tst));




disp(['Number of features: ' num2str(size(ufeats, 1))]);
disp('Train:');
disp(sprintf(' Specificity: %f +- %f', spec.tr));
disp(sprintf(' Sensitivity: %f +- %f\n', sens.tr));
disp('Test:');
disp(sprintf(' Specificity: %f +- %f', spec.tst));
disp(sprintf(' Sensitivity: %f +- %f\n', sens.tst));

experiment_name =  regexp(path, '(\d+)', 'match');
%save( char(strcat(path, '/', 'params_exp', experiment_name{:}, '_fixDur_', num2str(fixationDuration),'.mat')) , 'params' );

end
