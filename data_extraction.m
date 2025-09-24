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
        set_numbers = [2,10,13,15]; 
        sets = {set002, set010, set013, set015};
    case "Gaze"
        set_numbers = [3,8,14,18]; 
        sets = {set003, set008, set014, set018};
    case "Motor"
        set_numbers = [4,9,11,16];
        sets = {set004, set009, set011, set016};
    case "Gaze + Motor"
        set_numbers = [5,7,12,17]; 
        sets = {set005, set007, set012, set017};        
end 

get_start_idx = @(S) find(S.trial_num == 1, 1);
idx_start = arrayfun(@(i) get_start_idx(sets{i}), 1:numel(sets));

MaskedSpikeCounts = vertcat( ...
sets{1}.SpikeCount(idx_start(1):end, 1:5:end), ...
sets{2}.SpikeCount(idx_start(2):end, 1:5:end), ...
sets{3}.SpikeCount(idx_start(3):end, 1:5:end), ...
sets{4}.SpikeCount(idx_start(4):end, 1:5:end) ...
);

task_state_labels = vertcat( ...
sets{1}.TaskStateMasks.state_name(idx_start(1):end)', ...
sets{2}.TaskStateMasks.state_name(idx_start(2):end)', ...
sets{3}.TaskStateMasks.state_name(idx_start(3):end)', ...
sets{4}.TaskStateMasks.state_name(idx_start(4):end)' ...
);

target_coordinates = vertcat( ...
[sets{1}.TaskStateMasks.target(2, idx_start(1):end)'  sets{1}.TaskStateMasks.target(3, idx_start(1):end)'], ...
[sets{2}.TaskStateMasks.target(2, idx_start(2):end)'  sets{2}.TaskStateMasks.target(3, idx_start(2):end)'], ...
[sets{3}.TaskStateMasks.target(2, idx_start(3):end)'  sets{3}.TaskStateMasks.target(3, idx_start(3):end)'], ...
[sets{4}.TaskStateMasks.target(2, idx_start(4):end)'  sets{4}.TaskStateMasks.target(3, idx_start(4):end)'] ...
);

ends    = cellfun(@(S) S.trial_num(end), sets);
offsets = [0, cumsum(ends(1:end-1))];
trial_labels = vertcat( ...
sets{1}.trial_num(idx_start(1):end)' + offsets(1), ...
sets{2}.trial_num(idx_start(2):end)' + offsets(2), ...
sets{3}.trial_num(idx_start(3):end)' + offsets(3), ...
sets{4}.trial_num(idx_start(4):end)' + offsets(4) ...
);

if isequal( sort(string(sets{1}.TaskStateMasks.states)), ...
            sort(string(sets{2}.TaskStateMasks.states)), ...
            sort(string(sets{3}.TaskStateMasks.states)), ...
            sort(string(sets{4}.TaskStateMasks.states)) )            
    state_names = string(sets{1}.TaskStateMasks.states); 
end 

state_names(1) = 'Center';
empty_idx = cellfun(@isempty, task_state_labels);
task_state_labels(empty_idx) = {'Center'};
target_coordinates(empty_idx,1) = 0; 
target_coordinates(empty_idx,2) = 0; 

% Transform the "Remove" task state into "Reach" task state
for i = 1:numel(task_state_labels)   
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
trial_per_set = [40,32,32,32];  % da modificare
n_trials = sum(trial_per_set); 
% paradigm = {'Motor', 'Motor', 'Motor', 'Motor'};   % da modificare
paradigm = {'Free-gaze', 'Free-gaze', 'Free-gaze', 'Free-gaze'};
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
    data_by_trial = cell(numel(starts), 1); 
    mask_by_trial = cell(numel(starts), 1); 
    j = 1;
    for i = 1:numel(starts) 
        data_by_trial{j} = current_array(starts(i):ends(i), :);
        mask_by_trial{j,1}(:,1) = task_state_labels(starts(i):ends(i))';           
        mask_by_trial{j,1}(:,2) = num2cell(target_info(starts(i):ends(i),3)');
        j = j + 1; 
    end 

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
        target_id = vals(vals ~= 0);                            % da modificare                                                   
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
        len = [];        
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
clearvars -except dataset

save("free-gaze.mat", "dataset"); 
% load("free-gaze.mat")
% data_1 = dataset; 
% load("gaze.mat")
% data_2 = dataset; 
% load("motor.mat")
% data_3 = dataset;
% load("motor+gaze.mat");
% data_4 = dataset; 
% 
% dataset = [data_1; data_2; data_3; data_4]; 