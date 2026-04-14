clearvars
close all
clc

%% Baseline based on 60s quiet recording
[file, path] = uigetfile('*.mat', 'Seleziona un file MAT');
if isequal(file,0)
    error('Nessun file selezionato.');
else
    disp('Loading data...')
    load(fullfile(path, file));
end

counts = cellfun(@numel, baselineStruct.spikes.motor);
bas_m_medial = counts(1:96)./60;
bas_m_lateral = counts(97:192)./60;

baseline_const = [bas_m_medial, bas_m_lateral]; 

bin_size = 0.02; 
time_bins = 0:bin_size:60;
n_channels = 192;
for ch = 1:n_channels
    s = baselineStruct.spikes.motor{ch};
    tmp = histcounts(s, time_bins); 
    bin_baseline(: ,ch) = tmp/bin_size; 
end 

raw_des = length(dataset(1).Data(1).Resampled(1).Trial);
raw_act = size(bin_baseline,1);
nBlocks = floor(raw_act / raw_des);

for ch = 1:n_channels
    bin = bin_baseline(:,ch); 
    bin = bin(1 : nBlocks * raw_des);
    baseline(ch).mean =  mean(reshape(bin, raw_des, nBlocks), 2);
    baseline(ch).std = std(mean(reshape(bin, raw_des, nBlocks))); 
    baseline(ch).sem = std(mean(reshape(bin, raw_des, nBlocks)))/sqrt(length(mean(reshape(bin, raw_des, nBlocks)))); 
end 
clearvars -except dataset baseline baseline_const

%% Baseline as suggested by John
filename = {'../00_Data_extraction/free-gaze_BCI02_withtracker_exclUpdated.mat',... 
                   '../00_Data_extraction/motor_BCI02_withtracker_exclUpdated.mat',...
                   '../00_Data_extraction/controlled_BCI02_withtracker_exclUpdated.mat',...
                   '../00_Data_extraction/gaze_BCI02_withtracker_exclUpdated.mat'};

nCond = numel(filename);

% Single channel
n_sets = 6; 
sets_baseline = [1,2,3,4,5,6];
n_arrays = 2;
n_channels = 96;
n_trials = 32;  
bin_size = 0.02;

mean_baseline = nan(n_channels, n_arrays);
std_baseline  = nan(n_channels, n_arrays);
for array = 1:n_arrays
    for channel = 1:n_channels
        disp(array);
        disp(channel);
        trial_baselines = [];
        for d = 1:nCond
            load(filename{d});
            for set_ = 1:n_sets
                set = sets_baseline(set_);
                for trial = 1:n_trials
                     trial_spikes = data(set).Data(array).Interp(trial).Trial(:,channel);  
                     trial_baseline_counts = mean(trial_spikes); 
                     trial_baselines(end+1) = trial_baseline_counts;            
                end
            end 
        end
        baseline_mean_counts = mean(trial_baselines);
        baseline_std_counts  = std(trial_baselines);

        mean_baseline(channel, array) = baseline_mean_counts / bin_size; 
        std_baseline(channel, array)  = baseline_std_counts  / bin_size;
    end 
end 


%% Baseline as suggested by Charles - condition/targets
filename = {'../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02.mat',... 
           '../00_Data_extraction/BCI02_Session_0924/motor_BCI02.mat',...
           '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02.mat'};


% Single channel
sets_baseline = [1,2,3,4,5,6];
n_sets = length(sets_baseline);
n_arrays = 2;
n_channels = 96;
n_targets = 8; 
n_trials = 16;
bin_size = 0.02;
period = 0.1; 
bins_period = max(1, round(period/bin_size));


for d = 1:numel(filename)
    load(filename{d});

    n_ch = n_arrays*n_channels;
    n_tot_trials = n_sets*n_trials;

    % Baseline by trials
    bmean_by_trial = nan(n_tot_trials, n_ch);
    bstd_by_trial  = nan(n_tot_trials, n_ch);

    % Baseline by targets
    bmean_by_target = nan(n_targets, n_ch);
    bstd_by_target  = nan(n_targets, n_ch);

    phase_id = 1;
    for array = 1:n_arrays
        for channel = 1:n_channels
            ch_global = (array-1)*n_channels + channel;
            disp(array);
            disp(channel);
    
            trial_baselines_mean = nan(n_tot_trials, 1);
            trial_baselines_std = nan(n_tot_trials, 1);
            trial_targets_all = nan(n_tot_trials, 1);

            idx = 0;
            for set_ = 1:n_sets
                set = sets_baseline(set_);
                
                for trial = 1:n_trials
                    idx = idx + 1;
               
                    trial_spikes = data(set).Data(array).Interp(trial).Task_states{phase_id, 2}(end-bins_period+1:end,channel);  
                    trial_baselines_mean(idx) = mean(trial_spikes); 
                    trial_baselines_std(idx) = std(trial_spikes);

                    trial_targets_all(idx) = data(set).Data(array).Interp(trial).Target_ID;
                end 
            end
            bmean_by_trial(:,ch_global) = trial_baselines_mean./bin_size;
            bstd_by_trial(:,ch_global) = trial_baselines_std./bin_size; 

            for target = 1:n_targets
                mask = (trial_targets_all == target);

                if any(mask)
                    X = [];
                    idx = 0;
                    for set_ = 1:n_sets
                        set = sets_baseline(set_);

                        for trial = 1:n_trials
                            idx = idx + 1;
                            if ~mask(idx), continue; end
    
                            trial_spikes = data(set).Data(array).Interp(trial).Task_states{phase_id,2}(end-bins_period+1:end, channel);
                            X(:, end+1) = trial_spikes;
                            % X = [X; mean(trial_spikes,1)];
                        end
                    end
                    
                    mean_across_trials = mean(X, 2);                 
                    bmean_by_target(target, ch_global) = mean(mean_across_trials)./bin_size;  
                    bstd_by_target(target, ch_global) = std(mean_across_trials)./bin_size;
                    % bmean_by_target(target, ch_global) = mean(X)./bin_size;  
                    % bstd_by_target(target, ch_global) = std(X)./bin_size;
                end
            end
        end 
     end 

    mean_baseline.by_trials{1,d} = bmean_by_trial;
    std_baseline.by_trials{1,d} = bstd_by_trial;
    mean_baseline.by_targets{1,d} = bmean_by_target;
    std_baseline.by_targets{1,d} = bstd_by_target;
end 


%% Baseline globale
filenames = {
    '../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat', ...
    '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat', ...
    '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat', ...
    '../00_Data_extraction/BCI02_Session_0924/gaze_BCI02_exclUpdated.mat'
    };

n_sets     = 6;
n_arrays   = 2;
n_channels = 96;
n_trials   = 16;

bin_size    = 0.02;
period      = 0.1;
bins_period = max(1, round(period / bin_size));

phase_id = 1;
n_ch     = n_arrays * n_channels;

verbose = false;

% Carica i file una sola volta
all_data = cell(size(filenames));
for d = 1:numel(filenames)
    tmp = load(filenames{d});
    all_data{d} = tmp.data;
end

sum_baseline   = zeros(1, n_ch);
sumsq_baseline = zeros(1, n_ch);
count_baseline = zeros(1, n_ch);

for d = 1:numel(all_data)
    data = all_data{d};

    for set_idx = 1:n_sets
        for array = 1:n_arrays
            for trial = 1:n_trials

                % Salta trial esclusi
                if data(set_idx).Data(array).Interp(trial).Excluded ~= 0
                    if verbose
                        fprintf('Excluded trial: file=%d, set=%d, array=%d, trial=%d\n', d, set_idx, array, trial);
                    end
                    continue;
                end

                x = data(set_idx).Data(array).Interp(trial).Task_states{phase_id, 2};

                if isempty(x)
                    continue;
                end

                if size(x, 1) < bins_period
                    if verbose
                        fprintf('Too few bins: file=%d, set=%d, array=%d, trial=%d\n', d, set_idx, array, trial);
                    end
                    continue;
                end

                % Baseline: ultime bins_period righe, tutti i canali dell'array
                bs = x(end-bins_period+1:end, 1:n_channels);   % [bins_period x n_channels]

                % Indici globali dei canali dell'array corrente
                ch_idx = (array-1)*n_channels + (1:n_channels);

                % Statistiche con omitnan per colonna
                valid_mask = ~isnan(bs);
                count_baseline(ch_idx) = count_baseline(ch_idx) + sum(valid_mask, 1);
                sum_baseline(ch_idx)   = sum_baseline(ch_idx)   + sum(bs, 1, 'omitnan');
                sumsq_baseline(ch_idx) = sumsq_baseline(ch_idx) + sum(bs.^2, 1, 'omitnan');
            end
        end
    end
end

% Media
mean_baseline_common = sum_baseline ./ count_baseline;

% Deviazione standard campionaria
std_baseline_common = nan(1, n_ch);
valid_ch = count_baseline > 1;

std_baseline_common(valid_ch) = sqrt( ...
    (sumsq_baseline(valid_ch) - ...
    (sum_baseline(valid_ch).^2 ./ count_baseline(valid_ch))) ./ ...
    (count_baseline(valid_ch) - 1) );

% Conversione in firing rate
mean_baseline_common = mean_baseline_common / bin_size;
std_baseline_common  = std_baseline_common  / bin_size;



%% Baseline globale per condizione
filenames = {
    '../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat', ...
    '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat', ...
    '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat'
};

condition_names = {'free-gaze', 'motor', 'controlled','gaze'};

n_sets     = 6;
n_arrays   = 2;
n_channels = 96;
n_trials   = 16;

bin_size    = 0.02;
period      = 0.1;
bins_period = max(1, round(period / bin_size));

phase_id = 1;
n_ch     = n_arrays * n_channels;
n_cond   = numel(filenames);

verbose = false;

% Carica i file una sola volta
all_data = cell(size(filenames));
for d = 1:n_cond
    tmp = load(filenames{d});
    all_data{d} = tmp.data;
end

% Statistiche separate per condizione
sum_baseline   = zeros(n_cond, n_ch);
sumsq_baseline = zeros(n_cond, n_ch);
count_baseline = zeros(n_cond, n_ch);

for d = 1:n_cond
    data = all_data{d};

    for set_idx = 1:n_sets
        for array = 1:n_arrays
            for trial = 1:n_trials

                % Salta trial esclusi
                if data(set_idx).Data(array).Interp(trial).Excluded ~= 0
                    if verbose
                        fprintf('Excluded trial: cond=%d, set=%d, array=%d, trial=%d\n', ...
                            d, set_idx, array, trial);
                    end
                    continue;
                end

                x = data(set_idx).Data(array).Interp(trial).Task_states{phase_id, 2};

                if isempty(x)
                    continue;
                end

                if size(x, 1) < bins_period
                    if verbose
                        fprintf('Too few bins: cond=%d, set=%d, array=%d, trial=%d\n', ...
                            d, set_idx, array, trial);
                    end
                    continue;
                end

                % Baseline: ultime bins_period righe, tutti i canali dell'array
                bs = x(end-bins_period+1:end, 1:n_channels);   % [bins_period x n_channels]

                % Indici globali dei canali dell'array corrente
                ch_idx = (array-1)*n_channels + (1:n_channels);

                % Statistiche per questa condizione
                valid_mask = ~isnan(bs);
                count_baseline(d, ch_idx) = count_baseline(d, ch_idx) + sum(valid_mask, 1);
                sum_baseline(d, ch_idx)   = sum_baseline(d, ch_idx)   + sum(bs, 1, 'omitnan');
                sumsq_baseline(d, ch_idx) = sumsq_baseline(d, ch_idx) + sum(bs.^2, 1, 'omitnan');
            end
        end
    end
end

% Media baseline per condizione
mean_baseline_cond = sum_baseline ./ count_baseline;

% Deviazione standard campionaria per condizione
std_baseline_cond = nan(n_cond, n_ch);
valid_ch = count_baseline > 1;

std_baseline_cond(valid_ch) = sqrt( ...
    (sumsq_baseline(valid_ch) - ...
    (sum_baseline(valid_ch).^2 ./ count_baseline(valid_ch))) ./ ...
    (count_baseline(valid_ch) - 1) );

% Conversione in firing rate
mean_baseline_cond = mean_baseline_cond / bin_size;
std_baseline_cond  = std_baseline_cond  / bin_size;

% Facoltativo: output strutturato per condizione
baseline_stats = struct();
for d = 1:n_cond
    baseline_stats(d).condition = condition_names{d};
    baseline_stats(d).mean = mean_baseline_cond(d, :);
    baseline_stats(d).std  = std_baseline_cond(d, :);
    baseline_stats(d).count = count_baseline(d, :);
end