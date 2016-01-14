function [ sensivity,specificity] = optimal_bandpass(  eegTp, eegNTp, ufeats, A1_ch, A2_ch, fixationDuration, sRate, chans_for_clf )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    lowFq = 15;
    highFq = 50;
    band_length = 10;
    halfRate = sRate/2;
    for band = lowFq:(highFq - band_length)+1
        [b,a] = butter(5,[band,band+band_length-1]/halfRate, 'bandpass');
        eegTwp = filter(b,a,eegTp);
        eegNTwp = filter(b,a,eegNTp);
        [sens(band),spec(band)] = optimal_beta( eegTwp, eegNTwp, ufeats, A1_ch, A2_ch, fixationDuration, sRate, chans_for_clf );
        
    end;
    sens_tst_val=[];
    sens_tst_tau=[];
    sens_tst_dt=[];
    sens_tst_beta=[];
    
    spec_tst_val=[];
    spec_tst_tau=[];
    spec_tst_dt=[];
    spec_tst_beta=[];
    for i=lowFq:(highFq - band_length)+1
        sens_tst_val(i-lowFq+1) = sens(i).test.value(1);
        sens_tst_tau(i-lowFq+1) = sens(i).test.tau;
        sens_tst_dt(i-lowFq+1) = sens(i).test.dt;
        sens_tst_beta(i-lowFq+1) = sens(i).test.beta;
        
        spec_tst_val(i-lowFq+1) = spec(i).test.value(1);
        spec_tst_tau(i-lowFq+1) = spec(i).test.tau;
        spec_tst_dt(i-lowFq+1) = spec(i).test.dt;
        spec_tst_beta(i-lowFq+1) = spec(i).test.beta;
        
    end;
    subplot (4,1,1);
    plot(lowFq:(highFq - band_length)+1,sens_tst_val);
    title('sens value');
    subplot (4,1,2);
    plot(lowFq:(highFq - band_length)+1,sens_tst_tau);
    title('sens tau');
    subplot (4,1,3);
    plot(lowFq:(highFq - band_length)+1,sens_tst_dt);
    title('sens dt');
    subplot (4,1,4);
    plot(lowFq:(highFq - band_length)+1,sens_tst_beta);
    title('sens beta');
    
    figure(2),subplot (4,1,1);
    plot(lowFq:(highFq - band_length)+1,spec_tst_val);
    title('spec value');
    subplot (4,1,2);
    plot(lowFq:(highFq - band_length)+1,spec_tst_tau);
    title('spec tau');
    subplot (4,1,3);
    plot(lowFq:(highFq - band_length)+1,spec_tst_dt);
    title('spec dt');
    subplot (4,1,4);
    plot(lowFq:(highFq - band_length)+1,spec_tst_beta);
    title('spec beta');
end

