clear all
close all
clc

[file, path] = uigetfile('*.mat', 'Seleziona un file MAT');
if isequal(file,0)
    error('Nessun file selezionato.');
else
    disp('Loading data...')
    load(fullfile(path, file));
end


%% Select informative variables 

paradigm = input("Please enter the condition you want to analyse (i.e., 'Free-gaze', 'Gaze', 'Motor', 'Gaze + Motor'): ");
switch (paradigm) 
    case "Free-gaze"
        % set002, set010, set013, set015          
        set_numbers = [2,10,13,15]; 
        sets = {set002, set010, set013, set015};

        % idxs of time bins at the beginning of each set to be removed
        empty_idx = cellfun(@isempty, set002.TaskStateMasks.state_name');
        d = diff([empty_idx; 0]);
        idx_start_1 = find(d == -1);
        idx_start_1 = idx_start_1 - 1; 
        clear empty_idx

        empty_idx = cellfun(@isempty, set010.TaskStateMasks.state_name');
        d = diff([empty_idx; 0]);
        idx_start_2 = find(d == -1);
        idx_start_2 = idx_start_2 - 1; 
        clear empty_idx

        empty_idx = cellfun(@isempty, set013.TaskStateMasks.state_name');
        d = diff([empty_idx; 0]);
        idx_start_3 = find(d == -1);
        idx_start_3 = idx_start_3 - 1; 
        clear empty_idx

        empty_idx = cellfun(@isempty, set015.TaskStateMasks.state_name');
        d = diff([empty_idx; 0]);
        idx_start_4 = find(d == -1);
        idx_start_4 = idx_start_4 - 1; 
        clear empty_idx

        idx_start = [idx_start_1, idx_start_2, idx_start_3, idx_start_4];


        MaskedSpikeCounts = [set002.SpikeCount(idx_start(1):end,1:5:end); set010.SpikeCount(idx_start(2):end,1:5:end); set013.SpikeCount(idx_start(3):end,1:5:end); set015.SpikeCount(idx_start(4):end,1:5:end)]; 
        task_state_labels = [set002.TaskStateMasks.state_name(idx_start(1):end)'; set010.TaskStateMasks.state_name(idx_start(2):end)'; set013.TaskStateMasks.state_name(idx_start(3):end)'; set015.TaskStateMasks.state_name(idx_start(4):end)'];
        target_coordinates = [[set002.TaskStateMasks.target(2,idx_start(1):end)'; set010.TaskStateMasks.target(2,idx_start(2):end)'; set013.TaskStateMasks.target(2,idx_start(3):end)'; set015.TaskStateMasks.target(2,idx_start(4):end)'], ...
                            [set002.TaskStateMasks.target(3,idx_start(1):end)'; set010.TaskStateMasks.target(3,idx_start(2):end)'; set013.TaskStateMasks.target(3,idx_start(3):end)'; set015.TaskStateMasks.target(3,idx_start(4):end)']];
        trial_labels = [set002.trial_num(idx_start(1):end)'; set010.trial_num(idx_start(2):end)' + set002.trial_num(end); set013.trial_num(idx_start(3):end)' + set002.trial_num(end) + set010.trial_num(end); set015.trial_num(idx_start(4):end)' + set002.trial_num(end) + set010.trial_num(end) + set013.trial_num(end)]; 
        if isequal( sort(string(set002.TaskStateMasks.states)), ...
                    sort(string(set010.TaskStateMasks.states)), ...
                    sort(string(set013.TaskStateMasks.states)), ...
                    sort(string(set015.TaskStateMasks.states)) )            
            state_names = string(set002.TaskStateMasks.states); 
        end 
    case "Gaze"
        % set003, set008, set014, set018
        set_numbers = [3,8,14,18]; 

        % idxs of time bins at the beginning of each set to be removed
        empty_idx = cellfun(@isempty, set003.TaskStateMasks.state_name');
        d = diff([empty_idx; 0]);
        idx_start = find(d == -1);
        idx_start = idx_start - 1; 
        clear empty_idx

        empty_idx = cellfun(@isempty, set003.TaskStateMasks.state_name');


        MaskedSpikeCounts = [set003.SpikeCount(idx_start(1):end,1:5:end); set008.SpikeCount(idx_start(2):end,1:5:end); set014.SpikeCount(idx_start(3):end,1:5:end); set018.SpikeCount(idx_start(3):end,1:5:end)]; 
        task_state_labels = [set003.TaskStateMasks.state_name(idx_start(1):end)'; set008.TaskStateMasks.state_name(idx_start(2):end)'; set014.TaskStateMasks.state_name(idx_start(3):end)'; set018.TaskStateMasks.state_name(idx_start(4):end)'];
        target_coordinates = [[set003.TaskStateMasks.target(2,idx_start(1):end)'; set008.TaskStateMasks.target(2,idx_start(2):end)'; set014.TaskStateMasks.target(2,idx_start(3):end)'; set018.TaskStateMasks.target(2,idx_start(4):end)'], ...
                            [set003.TaskStateMasks.target(3,idx_start(1):end)'; set008.TaskStateMasks.target(3,idx_start(2):end)'; set014.TaskStateMasks.target(3,idx_start(3):end)'; set018.TaskStateMasks.target(3,idx_start(4):end)']];
        trial_labels = [set003.trial_num(idx_start(1):end)'; set008.trial_num(idx_start(2):end)' + set003.trial_num(end); set014.trial_num(idx_start(3):end)' + set003.trial_num(end) + set008.trial_num(end); set018.trial_num(idx_start(4):end)' + set003.trial_num(end) + set008.trial_num(end) + set014.trial_num(end)];     
        if isequal( sort(string(set003.TaskStateMasks.states)), ...
                    sort(string(set008.TaskStateMasks.states)), ...
                    sort(string(set014.TaskStateMasks.states)), ...
                    sort(string(set018.TaskStateMasks.states)) )            
            state_names = string(set003.TaskStateMasks.states); 
        end 
    case "Motor"
        % set004, set009, set011, set016
        set_numbers = [4,9,11,16];
        task_state_labels = [set004.TaskStateMasks.state_name'; set009.TaskStateMasks.state_name'; set011.TaskStateMasks.state_name'; set016.TaskStateMasks.state_name'];

        % idxs of time bins at the beginning of each set to be removed
        empty_idx = cellfun(@isempty, task_state_labels);
        d = diff([empty_idx; 0]);
        idx_start = find(d == -1);
        idx_start = idx_start - 1; 
        clear task_state_labels;

        MaskedSpikeCounts = [set004.SpikeCount(idx_start(1):end,1:5:end); set009.SpikeCount(idx_start(2):end,1:5:end); set011.SpikeCount(idx_start(3):end,1:5:end); set016.SpikeCount(idx_start(4):end,1:5:end)];
        task_state_labels = [set004.TaskStateMasks.state_name(idx_start(1):end)'; set009.TaskStateMasks.state_name(idx_start(2):end)'; set011.TaskStateMasks.state_name(idx_start(3):end)'; set016.TaskStateMasks.state_name(idx_start(4):end)'];
        target_coordinates = [[set004.TaskStateMasks.target(2,idx_start(1):end)'; set009.TaskStateMasks.target(2,idx_start(2):end)'; set011.TaskStateMasks.target(2,idx_start(3):end)'; set016.TaskStateMasks.target(2,idx_start(4):end)'], ...
                            [set004.TaskStateMasks.target(3,idx_start(1):end)'; set009.TaskStateMasks.target(3,idx_start(2):end)'; set011.TaskStateMasks.target(3,idx_start(3):end)'; set016.TaskStateMasks.target(3,idx_start(4):end)']];
        trial_labels = [set004.trial_num(idx_start(1):end)'; set009.trial_num(idx_start(2):end)' + set004.trial_num(end); set011.trial_num(idx_start(3):end)' + set009.trial_num(end) + set004.trial_num(end); set016.trial_num(idx_start(4):end)' + set011.trial_num(end) + set009.trial_num(end) + set004.trial_num(end)]; 
        if isequal( sort(string(set004.TaskStateMasks.states)), ...
                    sort(string(set009.TaskStateMasks.states)), ...
                    sort(string(set011.TaskStateMasks.states)), ...
                    sort(string(set016.TaskStateMasks.states)) )            
            state_names = string(set004.TaskStateMasks.states); 
        end 
    case "Gaze + Motor"
        % set005, set007, set012, set017
        set_numbers = [5,7,12,17]; 
        task_state_labels = [set005.TaskStateMasks.state_name'; set007.TaskStateMasks.state_name'; set012.TaskStateMasks.state_name'; set017.TaskStateMasks.state_name'];

        % idxs of time bins at the beginning of each set to be removed
        empty_idx = cellfun(@isempty, task_state_labels);
        d = diff([empty_idx; 0]);
        idx_start = find(d == -1);
        idx_start = idx_start - 1; 
        clear task_state_labels;

        MaskedSpikeCounts = [set005.SpikeCount(idx_start(1):end,1:5:end); set007.SpikeCount(idx_start(2):end,1:5:end); set012.SpikeCount(idx_start(3):end,1:5:end); set017.SpikeCount(idx_start(4):end,1:5:end)];
        task_state_labels = [set005.TaskStateMasks.state_name(idx_start(1):end)'; set007.TaskStateMasks.state_name(idx_start(2):end)'; set012.TaskStateMasks.state_name(idx_start(3):end)'; set017.TaskStateMasks.state_name(idx_start(4):end)'];        
        target_coordinates = [[set005.TaskStateMasks.target(2,idx_start(1):end)'; set007.TaskStateMasks.target(2,idx_start(2):end)'; set012.TaskStateMasks.target(2,idx_start(3):end)'; set017.TaskStateMasks.target(2,idx_start(4):end)'], ...
                            [set005.TaskStateMasks.target(3,idx_start(1):end)'; set007.TaskStateMasks.target(3,idx_start(2):end)'; set012.TaskStateMasks.target(3,idx_start(3):end)'; set017.TaskStateMasks.target(3,idx_start(4):end)']];
        trial_labels = [set005.trial_num(idx_start(1):end)'; set007.trial_num(idx_start(2):end)' + set005.trial_num(end); set012.trial_num(idx_start(3):end)' + set007.trial_num(end) + set005.trial_num(end); set017.trial_num(idx_start(4):end)' set012.trial_num(end) + set007.trial_num(end) + set005.trial_num(end)]; 
        if isequal( sort(string(set005.TaskStateMasks.states)), ...
                    sort(string(set007.TaskStateMasks.states)), ...
                    sort(string(set012.TaskStateMasks.states)), ...
                    sort(string(set017.TaskStateMasks.states)) )            
            state_names = string(set005.TaskStateMasks.states); 
        end 
end 

state_names(1) = 'Center';
empty_idx = cellfun(@isempty, task_state_labels);
task_state_labels(empty_idx) = {'Center'};
target_coordinates(empty_idx,1) = 0; 
target_coordinates(empty_idx,2) = 0; 

% Transform the "Remove" task state into "Reach" task state
for i = 2:numel(task_state_labels)   
    if strcmp(task_state_labels{i}, 'Remove') | strcmp(task_state_labels{i}, 'Remove1') | strcmp(task_state_labels{i}, 'Remove2')
        task_state_labels{i} = 'Reach';
    end
end

% Removing "Remove" from the possible state names
if any(state_names == "Remove" | state_names == "Remove1" | state_names == "Remove2") 
    state_names(state_names == "Remove" | state_names == "Remove1" | state_names == "Remove2") = [];
end

% Assigning labels to targets
target_info = target_coordinates; 
target_info(:,3) = NaN;
refXY = [  0      0
           0      0.2
           0.14   0.14
           0.2    0
           0.14  -0.14
           0     -0.2
          -0.14  -0.14
          -0.2    0
          -0.14   0.14 ];
codes = (0:8).';
tol = 1e-6;
[tf, loc] = ismembertol(target_info(:,1:2), refXY, tol, 'ByRows', true);
target_info(tf,3) = codes(loc(tf));

%% Distinction between arrays
% MaskedSpikeCounts = Data.SpikeCount(:,1:5:end);
rowsWithNaN = any(isnan(MaskedSpikeCounts), 2);
rowsToRemove = find(rowsWithNaN);
MaskedSpikeCounts(rowsToRemove, :) = [];
task_state_labels(rowsToRemove) = [];
target_coordinates(rowsToRemove,:) = [];
trial_labels(rowsToRemove) = []; 

idx.medial_sens   = 65:96;                   % Sensory medial (32ch)
idx.lateral_sens  = 193:224;                 % Sensory lateral (32ch)
idx.medial_motor  = [1:64,   97:128];        % Motor medial (96ch)
idx.lateral_motor = [129:192, 225:256];      % Motor lateral (96ch)

arrays.motor_m = MaskedSpikeCounts(:, idx.medial_motor);
arrays.motor_l = MaskedSpikeCounts(:, idx.lateral_motor);
arrays.sens_m  = MaskedSpikeCounts(:, idx.medial_sens);
arrays.sens_l  = MaskedSpikeCounts(:, idx.lateral_sens);

%% Data structure definition
trial_per_set = [40,32,32,32]; % da modificare
n_trials = sum(trial_per_set); 
paradigm = {'Gaze', 'Gaze', 'Gaze', 'Gaze'}; % da modificare
struct_l1 = {'Set', 'Paradigm', 'Data'};
dataset_l1 = cell2struct(repmat({[]}, 1, numel(struct_l1)), struct_l1, 2);
n = numel(set_numbers);
dataset_l1 = repmat(dataset_l1, n, 1);   


%% Dataset based on trials, task states, target IDs 
starts = [1; find(diff(trial_labels(:)) ~= 0) + 1];
ends   = [starts(2:end) - 1; size(MaskedSpikeCounts,1)];

array_names = fieldnames(arrays);
for array = 1:numel(array_names)
    disp("Array: " + array_names{array})
    current_array = arrays.(array_names{array});

    % Trial split
    j = 1;
    data_by_trial = cell(numel(starts), 1); 
    mask_by_trial = cell(numel(starts), 2); 

    for i = 1:numel(starts) 
        data_by_trial_tmp = current_array(starts(i):ends(i), :);
        data_by_trial{j} = current_array(starts(i):ends(i), :);
        mask_by_trial{j,1}(:,1) = task_state_labels(starts(i):ends(i))';           
        mask_by_trial{j,1}(:,2) = num2cell(target_info(starts(i):ends(i),3)');
        j = j + 1; 
    end
    data_by_trial(j:end) = [];
    mask_by_trial(j:end) = [];
    mask_by_trial = mask_by_trial'; 

    % Task phases split
    data_by_task_state = cell(length(mask_by_trial),1);
    for j = 1:length(mask_by_trial)
        data_by_task_state_tmp = {};
        for i = 1:length(state_names)
            idx_tmp = find(strcmp(string(mask_by_trial{j,1}(:,1)), state_names(i)) == 1);
            
            if state_names(i) == "Center"                                                       % da modificare 
                d = diff(idx_tmp);
                breaks = [0 find(d > 1) length(idx_tmp)];
                blocks = cell(length(breaks)-1,1);
                for k = 1:length(breaks)-1
                    blocks{k} = idx_tmp(breaks(k)+1 : breaks(k+1));
                end
                if i == 1
                    idx_tmp = blocks{1,1};
                else
                    idx_tmp = blocks{2,1};
                end 
            end 

            if ~isempty(idx_tmp)
                data_by_task_state_tmp = [data_by_task_state_tmp; {state_names(1,i), data_by_trial{j,1}(idx_tmp, :)}];
            end        
        end  
        data_by_task_state{j} = data_by_task_state_tmp; 
        vals = cell2mat(mask_by_trial{j,1}(:,2));
        target_id = vals(vals ~= 0);      % da modificare                                                   
        has_two = numel(unique(target_id)) >= 2;                                                 
        if has_two  
            warning('target_id contains at least two different numbers.');
            disp(j);
        elseif isempty(target_id)
            warning('target_id is empty after filtering.');
            disp(j);
        end 

        if isempty(target_id)                                                                   
            target_ids_by_trial(j) = 0; 
        else 
            target_ids_by_trial(j) = target_id(1);
        end 
    end 

    % Data interpolation 
    nPhases = numel(state_names);
    nTrials = numel(data_by_task_state);
    data_by_task_state_int = data_by_task_state;

    for phase = 1:nPhases
        len = zeros(nTrials,1);        
        for trial = 1:nTrials
            len = [len; size(cell2mat(data_by_task_state{trial,1}(phase,2)),1)];
        end
        [maxLen, idx] = max(len); 

        for trial = 1:nTrials
            data_by_task_state_int{trial,1}(phase,2) = {interp_data(cell2mat(data_by_task_state{trial,1}(phase,2)), maxLen)};
        end
    end

    data_by_trial_int = data_by_trial; 
    for trial = 1:nTrials
        data_by_trial_int_tmp = []; 
        for phase = 1:nPhases
            data_by_trial_int_tmp = [data_by_trial_int_tmp; cell2mat(data_by_task_state_int{trial,1}(phase,2))];
        end 
        data_by_trial_int{trial, 1} = data_by_trial_int_tmp; 
    end 
    
    % Filling data structure 
    start_set = 1; 
    for n_set = 1:length(set_numbers)

        struct_l2 = {'Trial','Task_states','Target_ID'};
        dataset_l2 = cell2struct(repmat({[]}, 1, numel(struct_l2)), struct_l2, 2);
        dataset_l2 = repmat(dataset_l2, trial_per_set(n_set), 1);

        k = start_set;
        for n_trial = 1:trial_per_set(n_set)
            dataset_l2(n_trial).Trial       = data_by_trial{k};
            dataset_l2(n_trial).Task_states = data_by_task_state{k}; 
            dataset_l2(n_trial).Target_ID   = target_ids_by_trial(k); 
            k = k + 1; 
        end 

        original_l2 = dataset_l2; 
 
        dataset_l2 = cell2struct(repmat({[]}, 1, numel(struct_l2)), struct_l2, 2);
        dataset_l2 = repmat(dataset_l2, trial_per_set(n_set), 1);        
        
        k = start_set;
        for n_trial = 1:trial_per_set(n_set)
            dataset_l2(n_trial).Trial       = data_by_trial_int{k};
            dataset_l2(n_trial).Task_states = data_by_task_state_int{k}; 
            dataset_l2(n_trial).Target_ID   = target_ids_by_trial(k); 
            k = k + 1;
        end 

        resampled_l2 = dataset_l2;

        start_set = k;

        per_set_data{n_set}(array).Array      = array_names{array};
        per_set_data{n_set}(array).Original   = original_l2;
        per_set_data{n_set}(array).Resampled  = dataset_l2;
    end 
end

for n_set = 1:numel(set_numbers)
    dataset_l1(n_set).Set      = string(set_numbers(n_set));
    dataset_l1(n_set).Paradigm = paradigm{n_set};
    dataset_l1(n_set).Data     = per_set_data{n_set};  
end

dataset = dataset_l1;
clearvars -except Data dataset