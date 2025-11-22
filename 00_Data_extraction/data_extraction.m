clear all
% close all
clc

[file, path] = uigetfile('*.mat', 'Seleziona un file MAT');
if isequal(file,0)
    error('Nessun file selezionato.');
else
    disp('Loading data...')
    load(fullfile(path, file));
end


%% Reading data from the dataset
paradigm = input('Please enter the condition you want to analyse (i.e., ''Free-gaze'', ''Gaze'', ''Motor'', ''Gaze + Motor'', ''New dataset''): ', 's');
paradigm = string(paradigm);
switch (paradigm) 
    case "Free-gaze"
        % BCI03 1 dataset
        % set_numbers = [2,10,13,15]; 
        % sets = {set002, set010, set013, set015};
        % BCI03 2 dataset
        % set_numbers = [4, 7, 9, 16, 20, 22]; 
        % sets = {set04, set07, set09, set16, set20, set22};
        % BCI02 
        set_numbers = [4, 9, 12, 14, 21, 25]; 
        sets = {set04, set09, set12, set14, set21, set25};
    case "Gaze"
        % BCI03 1 dataset
        % set_numbers = [3,8,14,18]; 
        % sets = {set003, set008, set014, set018};
        % BCI03 2 dataset
        % set_numbers = [2, 8, 12, 13, 18, 23]; 
        % sets = {set02, set08, set12, set13, set18, set23};
        % BCI02
        set_numbers = [1, 6, 13, 17, 18, 23]; 
        sets = {set01, set06, set13, set17, set18, set23};
    case "Motor"
        % BCI03 1 dataset        
        % set_numbers = [4,9,11,16];
        % sets = {set004, set009, set011, set016};
        % BCI03 2 dataset        
        % set_numbers = [3, 5, 10, 15, 17, 24];
        % sets = {set03, set05, set10, set15, set17, set24};       
        % BCI02
        set_numbers = [2, 7, 10, 15, 20, 22]; 
        sets = {set02, set07, set10, set15, set20, set22};
    case "Gaze + Motor"
        % BCI03 1 dataset                
        % set_numbers = [5,7,12,17]; 
        % sets = {set005, set007, set012, set017}; 
        % BCI03 2 dataset                
        % set_numbers = [1, 6, 11, 14, 19, 21];
        % sets = {set01, set06, set11, set14, set19, set21};    
        % BCI02
        set_numbers = [3, 5, 11, 16, 19, 24]; 
        sets = {set03, set05, set11, set16, set19, set24};
    case "New dataset"
        set_numbers = [1, 2, 3, 4];
        sets = {set01, set02, set03, set04}; 
end 

get_start_idx = @(S) find(S.trial_num == 1, 1);
idx_start = arrayfun(@(i) get_start_idx(sets{i}), 1:numel(sets));

ends    = cellfun(@(S) S.trial_num(end), sets);
offsets = [0, cumsum(ends(1:end-1))];

MaskedSpikeCounts_tmp = []; 
task_state_labels_tmp = {}; 
target_coordinates_tmp = []; 
trial_labels_tmp = []; 
for i = 1:length(set_numbers)
    MaskedSpikeCounts = [MaskedSpikeCounts_tmp; sets{i}.SpikeCount(idx_start(i):end, 1:5:end)];
    MaskedSpikeCounts_tmp = MaskedSpikeCounts; 
    
    task_state_labels = [task_state_labels_tmp; sets{i}.TaskStateMasks.state_name(idx_start(i):end)'];
    task_state_labels_tmp = task_state_labels; 

    target_coordinates = [target_coordinates_tmp; sets{i}.TaskStateMasks.target(2, idx_start(i):end)'  sets{i}.TaskStateMasks.target(3, idx_start(i):end)']; 
    target_coordinates_tmp = target_coordinates; 

    trial_labels = [trial_labels_tmp; sets{i}.trial_num(idx_start(i):end)' + offsets(i)]; 
    trial_labels_tmp = trial_labels; 

end 

states = cellfun(@(s) sort(string(s.TaskStateMasks.states(:))), sets, 'UniformOutput', false);
ref = states{1};
ok  = all(cellfun(@(v) isequal(v, ref), states));

if ok
    state_names = string(sets{1}.TaskStateMasks.states);
else
    disp("Problem!")
end

state_names(1) = 'Center';
empty_idx = cellfun(@isempty, task_state_labels);
task_state_labels(empty_idx) = {'Center'};
target_coordinates(empty_idx,1) = 0; 
target_coordinates(empty_idx,2) = 0; 

% Transform the "Remove" task state into "Reach" task state
for i = 1:numel(task_state_labels)   
    if strcmp(task_state_labels{i}, 'Remove') | strcmp(task_state_labels{i}, 'Remove1') | strcmp(task_state_labels{i}, 'Remove2')
       task_state_labels{i} = 'Center';
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
rowsWithNaN = any(isnan(MaskedSpikeCounts), 2);
rowsToRemove = find(rowsWithNaN);
MaskedSpikeCounts(rowsToRemove, :) = [];
task_state_labels(rowsToRemove) = [];
target_info(rowsToRemove,:) = [];
trial_labels(rowsToRemove) = []; 

idx.medial_sens   = 65:96;                   % Sensory medial (32ch)
idx.lateral_sens  = 193:224;                 % Sensory lateral (32ch)
idx.medial_motor  = [1:64,   97:128];        % Motor medial (96ch)
idx.lateral_motor = [129:192, 225:256];      % Motor lateral (96ch)

arrays.motor_m = MaskedSpikeCounts(:, idx.medial_motor);
arrays.motor_l = MaskedSpikeCounts(:, idx.lateral_motor);
arrays.sens_m  = MaskedSpikeCounts(:, idx.medial_sens);
arrays.sens_l  = MaskedSpikeCounts(:, idx.lateral_sens);

%% Initial data structure definition
trial_per_set = ends; 
n_trials = sum(trial_per_set); 
n_sets = numel(set_numbers);
paradigm = repmat(paradigm, 1, n_sets);
struct_l1 = {'Set', 'Paradigm', 'Data'};
dataset_l1 = cell2struct(repmat({[]}, 1, numel(struct_l1)), struct_l1, 2);
dataset_l1 = repmat(dataset_l1, n_sets, 1);   


%% Dataset based on trials, task states, target IDs 
starts = [1; find(diff(trial_labels(:)) ~= 0) + 1];
ends   = [starts(2:end) - 1; size(MaskedSpikeCounts,1)];

array_names = fieldnames(arrays);
for array = 1:numel(array_names)
    disp("Array: " + array_names{array})
    current_array = arrays.(array_names{array});

    % Trial split
    % output vector: data_by_trial 
    data_by_trial = cell(numel(starts), 1); % Contains the spike count matrix devided per trial
    mask_by_trial = cell(numel(starts), 1); % First column: task phase; second column: target id
    j = 1;
    for i = 1:numel(starts) 
        data_by_trial{j} = current_array(starts(i):ends(i), :);
        mask_by_trial{j,1}(:,1) = task_state_labels(starts(i):ends(i))';           
        mask_by_trial{j,1}(:,2) = num2cell(target_info(starts(i):ends(i),3)');
        j = j + 1; 
    end 

    % Task phases split
    % output vectors: 
    % data_by_task_state
    % target_ids_by_trial
    n_phases = numel(state_names);
    data_by_task_state = cell(length(mask_by_trial),1);
    for trial = 1:n_trials
        data_by_task_state_tmp = {};

        % After this cycle, we obtain a vector in which each row corresponds to a 
        % trial, where the spike count matrix is divided according to the task phase.
        for phase = 1:n_phases
            idx_tmp = find(strcmp(string(mask_by_trial{trial,1}(:,1)), state_names(phase)) == 1);
            
            % We have two center phases, at the beginning and end of the
            % task. For the first occurance of the center label in state_names
            % vector, we consider the first block of ids; for the last element 
            % of the state_names vector, we consider the second block.
            if state_names(phase) == "Center"                                                      
                d = diff(idx_tmp);
                breaks = [0 find(d > 1) length(idx_tmp)];
                blocks = cell(length(breaks)-1,1);
                for k = 1:length(breaks)-1
                    blocks{k} = idx_tmp(breaks(k)+1 : breaks(k+1));
                end
                if phase == 1 
                    idx_tmp = blocks{1,1};
                else
                    idx_tmp = blocks{2,1};
                end 
            end 
            data_by_task_state_tmp = [data_by_task_state_tmp; {state_names(1,phase), data_by_trial{trial,1}(idx_tmp, :)}];       
        end  
        data_by_task_state{trial} = data_by_task_state_tmp; 
        
        % Here we retrieve the target ID associated with each trial. According to our 
        % labeling, zero represents the central point. In a trial, we have 0 at the 
        % beginning, then the target ID, and 0 again at the end. 
        % We search for numbers different from 0.
        vals = cell2mat(mask_by_trial{trial,1}(:,2));
        target_id = vals(vals ~= 0);                                                                              
        has_two = numel(unique(target_id)) >= 2;                                                 
        if has_two  
            warning('Two targets associated to one trial!');
            disp('Trial ID:')
            disp(trial);
        elseif isempty(target_id)
            warning('No targets associated to the trial!');
            disp('Trial ID:')
            disp(trial);
        end 
        target_ids_by_trial(trial) = unique(target_id);
    end 

    % Data interpolation 
    % Interpolation of spike count matrix by column in order to have the 
    % different phases of the task with the same number of samples
    % output vectors: 
    % data_by_trial_int
    % data_by_task_state_int
    data_by_task_state_int = cell(size(data_by_task_state));
    
    for t = 1:numel(data_by_task_state)
        n_rows = size(data_by_task_state{t},1);
        data_by_task_state_int{t} = cell(n_rows,2);
        data_by_task_state_int{t}(:,1) = data_by_task_state{t}(:,1);  % riporto i nomi delle fasi
    end

    for phase = 1:n_phases
        len = [];        
        for trial = 1:n_trials
            len = [len; size(cell2mat(data_by_task_state{trial,1}(phase,2)),1)];
        end
        [maxLen, idx] = max(len); 
        fprintf("Phase %s. Max length: %d. Average length: %d.\n", state_names(phase), maxLen, round(mean(len)))
        if maxLen > round(mean(len)) + 5
            fprintf("Phase %s. Outlier: %d, with an average length of %d. Trial %d, Set %d\n", state_names(phase), maxLen, round(mean(len)), mod(idx-1, 32) + 1, ceil(idx / 32));
        end

        for trial = 1:n_trials
            data_by_task_state_int{trial,1}(phase,2) = {interp_data(cell2mat(data_by_task_state{trial,1}(phase,2)), maxLen)};
        end
    end

    data_by_trial_int = cell(size(data_by_trial));
    for trial = 1:n_trials
        data_by_trial_int_tmp = []; 
        for phase = 1:n_phases
            data_by_trial_int_tmp = [data_by_trial_int_tmp; cell2mat(data_by_task_state_int{trial,1}(phase,2))];
        end 
        data_by_trial_int{trial} = data_by_trial_int_tmp; 
    end 
    
    % Filling data structure 
    start_set = 1; 
    for set = 1:n_sets

        struct_l2 = {'Trial','Task_states','Target_ID'};
        dataset_l2 = cell2struct(repmat({[]}, 1, numel(struct_l2)), struct_l2, 2);
        dataset_l2 = repmat(dataset_l2, trial_per_set(set), 1);

        k = start_set;
        for n_trial = 1:trial_per_set(set)
            dataset_l2(n_trial).Trial       = data_by_trial{k};
            dataset_l2(n_trial).Task_states = data_by_task_state{k}; 
            dataset_l2(n_trial).Target_ID   = target_ids_by_trial(k); 
            k = k + 1; 
        end 

        original_l2 = dataset_l2; 
 
        dataset_l2 = cell2struct(repmat({[]}, 1, numel(struct_l2)), struct_l2, 2);
        dataset_l2 = repmat(dataset_l2, trial_per_set(set), 1);        
        
        k = start_set;
        for n_trial = 1:trial_per_set(set)
            dataset_l2(n_trial).Trial       = data_by_trial_int{k};
            dataset_l2(n_trial).Task_states = data_by_task_state_int{k}; 
            dataset_l2(n_trial).Target_ID   = target_ids_by_trial(k); 
            k = k + 1;
        end 

        resampled_l2 = dataset_l2;

        start_set = k;

        per_set_data{set}(array).Array      = array_names{array};
        per_set_data{set}(array).Original   = original_l2;
        per_set_data{set}(array).Interp  = dataset_l2;
    end 
end

for set = 1:n_sets
    dataset_l1(set).Set      = string(set_numbers(set));
    dataset_l1(set).Paradigm = paradigm(set);
    dataset_l1(set).Data     = per_set_data{set};  
end

data = dataset_l1;
clearvars -except data

save("free-gaze_BCI02.mat", "data"); 