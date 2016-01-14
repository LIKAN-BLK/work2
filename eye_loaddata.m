function [eegT, eegNT, parameters, fixationDuration, sRate, path, epoch_size] = eye_loaddata(path, fixation_threshold, epoch_size, left_border)
%
% make target and non-target EEG trials from data in files
%

% ---------- load data -------------------
json = fileread([path '/meta.json']);
[str ~] = parse_json(json);

epoch_size = epoch_size - left_border;
eegT = [];
eegNT = [];
for i = 1:length(str{1}.valid_files)
    edf_fname = str{1}.valid_files{i}.name_edf;
    ascii_fname = ['eventsLatencies' edf_fname '.ascii'];
    dat_fname = str{1}.valid_files{i}.name_dat;

    [signal, states, parameters] = load_bcidat([path '/' dat_fname]);
    tmp = find(signal(:,end) ~= 0);
    signal = signal(tmp(3)+1:end, :);
    
    sRate = parameters.SamplingRate.NumericValue;
    
    tmp = fileread([ path '/' [edf_fname '.asc'] ] );
    [T, M] = regexp(tmp, 'fixationDuration":([0-9]+)','tokens', 'match');
    fixationDuration = str2num(T{1}{1});
    
    if (fixation_threshold ~= fixationDuration)
        continue
    end

    fid = fopen([path '/' ascii_fname]);
    tmp = fgetl(fid); % read header
    B = textscan(fid, '%f %s'); % B{1} - times, B{2} - labels
    fclose(fid);

    [eegTi, eegNTi] = loaddata1(signal, B, epoch_size, sRate, left_border);    
    eegT = cat(3, eegT, eegTi);
    eegNT = cat(3, eegNT, eegNTi);
end

fixationDuration = fixation_threshold;
bline_t0 = 200;
bline_t1 = 300;
eegTp = eegT - repmat(mean(eegT(bline_t0*sRate/1000:bline_t1*sRate/1000,:,:),1), size(eegT,1), 1, 1);
eegNTp = eegNT - repmat(mean(eegNT(bline_t0*sRate/1000:bline_t1*sRate/1000,:,:),1), size(eegNT,1), 1, 1);
end

function [eegT, eegNT] = loaddata1(signal, B, epoch_size, sRate, left_border)
% ----------------------------------------

msgbuttonPressed_idx = find(strcmp(B{2}, 'msgbuttonPressed'));
msgbuttonPressed_times = B{1}(msgbuttonPressed_idx);
msgbuttonPressed_t = round(msgbuttonPressed_times * sRate);

msgballChosen_idx = find(strcmp(B{2}, 'msgballChosen'));
msgballChosen_times = B{1}(msgballChosen_idx);
msgballChosen_t = round(msgballChosen_times * sRate);

msgBallMoved_idx = find(strcmp(B{2}, 'msgBallMoved'));
msgBallMoved_times = B{1}(msgBallMoved_idx);
msgBallMoved_t = round(msgBallMoved_times * sRate);

msgClickedInBlockMode_idx = find(strcmp(B{2}, 'msgClickedInBlockMode'));
msgClickedInBlockMode_times = B{1}(msgClickedInBlockMode_idx);
msgClickedInBlockMode_t = round(msgClickedInBlockMode_times * sRate);

msgBallClickedInBlockedMode_idx = find(strcmp(B{2}, 'msgBallClickedInBlockedMode'));
msgBallClickedInBlockedMode_times = B{1}(msgBallClickedInBlockedMode_idx);
msgBallClickedInBlockedMode_t = round(msgBallClickedInBlockedMode_times * sRate);

msgBoardClickedInBlockedMode_idx = find(strcmp(B{2}, 'msgBoardClickedInBlockedMode'));
msgBoardClickedInBlockedMode_times = B{1}(msgBoardClickedInBlockedMode_idx);
msgBoardClickedInBlockedMode_t = round(msgBoardClickedInBlockedMode_times * sRate);

fixationDuration_t = round((epoch_size/1000) * sRate);

eventsT_t = [msgbuttonPressed_t; msgballChosen_t; msgBallMoved_t];
eventsNT_t = [msgClickedInBlockMode_t; msgBallClickedInBlockedMode_t; msgBoardClickedInBlockedMode_t];

% ---------- make EEG trials ------------------------------
nChannels = size(signal, 2);
%signal = signal - repmat(mean(signal, 2), 1, nChannels);
left_border = left_border/1000*sRate;
eegT = zeros(fixationDuration_t, nChannels, length(eventsT_t));
for i = 1:length(eventsT_t)
    eegT(:, :, i) = signal(eventsT_t(i)+left_border+1:eventsT_t(i)+fixationDuration_t+left_border, :);    
end

eegNT = zeros(fixationDuration_t, nChannels, length(eventsNT_t));
for i = 1:length(eventsNT_t)
    eegNT(:, :, i) = signal(eventsNT_t(i)+left_border+1:eventsNT_t(i)+fixationDuration_t+left_border, :);    
end
i;
end