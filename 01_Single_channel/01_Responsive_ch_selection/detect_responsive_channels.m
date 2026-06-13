%% Responsive channel detection

% This script identifies task-responsive channels by analyzing firing-rate
% modulation around the selected event (e.g., movement or gaze-shift go cue).
% Neural activity is compared against a baseline, and channels are classified
% as responsive when they exhibit a sustained excitatory or inhibitory
% modulation exceeding predefined amplitude and duration thresholds.
% Results are stored per target and summarized across arrays and channels.

clearvars -except mean_baseline_common std_baseline_common mean_baseline std_baseline 
close all
clc

%% Parameters definition
n_sets = 6;
n_arrays = 2;
n_channels = 96;
ch_start = 1;
ch_end = 96;
target_des = 1:8;
bin_size = 0.02;
w = 15; % smoothing parameter

% To be changed based on the experimental condition:
% filename = '../../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat'; 
% filename = '../../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat';
filename = '../../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat';
array_names = ["medial", "lateral"];

% Modulation search window:
% 'planning'  = Go cue - 0.5 s --> Go cue
% 'reaching'  = Go cue         --> Go cue + 3.0 s
% 'full'      = Go cue - 0.5 s --> Go cue + 3.0 s
% I search for modulation around the movement go cue for the cursor movement
% conditions (free-gaze, gaze-on-center, gaze-on-target); additionally,
% for gaze-on-target and gaze-only conditions, I can also search for
% modulation around the gaze-shift go cue. The intervals above in the comment
% refer to the first case.
response_window_mode = 'full';

% absolute threshold relative to baseline
min_delta_from_baseline = 12;   % Hz

% minimum modulation duration
min_duration_s = 0.06;          % 60 ms, 3 bins
min_duration_bins = ceil(min_duration_s / bin_size);

% baseline mode:
% 'global'          -> uses mean_baseline_common
% 'target_specific' -> uses mean_baseline.by_targets(target, ch_global)
% 'prefer_target'   -> uses the target-specific baseline if available,
%                      otherwise falls back to the global baseline
% the global baseline takes into account the first 100 ms of all conditions
% (including gaze-only).
baseline_mode = 'global'; 


%% Check baseline
has_global_baseline = exist('mean_baseline_common', 'var') && ~isempty(mean_baseline_common);
has_target_baseline = exist('mean_baseline', 'var') && ~isempty(mean_baseline.by_targets);

if ~has_global_baseline && ~has_target_baseline
    error('No baseline available: mean_baseline_common and/or mean_baseline_by_targets are required.');
end

%% Data loading and searching window definition
S = load(filename);
data = S.data;

TS = data(1).Data(1).Interp(1).Task_states;
state_names = string(TS(:,1));
state_dur_s = cellfun(@(x) size(x,1) * bin_size, TS(:,2));

state_onset_s = [0; cumsum(state_dur_s(1:end-1))];
get_onset = @(name) state_onset_s(find(strcmpi(state_names, name), 1, 'first'));

% To be changed based on the experimental condition:
% task_name = "Free-gaze"; 
% lines_t   = [get_onset("pres12"); get_onset("reach")];
% lines_lab = ["Target cue"; "Go cue"];
%------
% task_name = "Gaze-on-center"; 
% lines_t   = [get_onset("pres12"); get_onset("reach")];
% lines_lab = ["Target cue"; "Go cue"];
%------
task_name = "Gaze-on-target"; 
lines_t   = [get_onset("Gaze"); get_onset("reach")];
lines_lab = ["Target cue"; "Go cue"];

idx_go = find(strcmp(lines_lab, "Go cue"), 1);
if isempty(idx_go)
    error('Go cue not found in the file''s Task_states.');
end
go_time = lines_t(idx_go);

% window definition
switch lower(response_window_mode)
    case 'planning'
        win_start = go_time - 0.5;
        win_end   = go_time;
    case 'reaching'
        win_start = go_time;
        win_end   = go_time + 3.0;
    case 'full'
        win_start = go_time - 0.5;
        win_end   = go_time + 3.0;
    otherwise
        error('Unrecognized response_window_mode.');
end

%% Output struct definition
responsive_channels = struct();
results = table();
responsive_matrix = false(numel(target_des), n_arrays*n_channels);

%% Target / Array / Channel loop
for t_id = 1:numel(target_des)
    target_id = target_des(t_id);

    % Initialization
    responsive_channels(t_id).target_id = target_id;
    responsive_channels(t_id).channels_global = [];
    responsive_channels(t_id).channels_local = [];
    responsive_channels(t_id).array_name = strings(0,1);
    responsive_channels(t_id).response_type = strings(0,1);
    responsive_channels(t_id).baseline_type = strings(0,1);

    fprintf('\n==============================\n');
    fprintf('TARGET %d\n', target_id);
    fprintf('================================\n');

    for array = 1:n_arrays
        for channel = ch_start:ch_end

            ch_global = (array - 1) * n_channels + channel;

            % --------------------------------------------------------
            % Baseline selection
            % --------------------------------------------------------
            baseline_val = NaN;
            baseline_type_used = "";

            % Case global
            baseline_global_val = NaN;
            if has_global_baseline && ch_global <= numel(mean_baseline_common)
                baseline_global_val = mean_baseline_common(ch_global);
            end

            % Case target-specific
            baseline_target_val = NaN;
            if has_target_baseline
                if target_id <= size(mean_baseline.by_targets{1,1},1) && ch_global <= size(mean_baseline.by_targets{1,1},2)
                    baseline_target_val = mean_baseline.by_targets{1,1}(target_id, ch_global);
                    if baseline_target_val == 0
                        baseline_target_val = NaN;
                    end
                end
            end

            switch lower(baseline_mode)
                case 'global'
                    baseline_val = baseline_global_val;
                    baseline_type_used = "global";

                case 'target_specific'
                    baseline_val = baseline_target_val;
                    baseline_type_used = "target_specific";

                case 'prefer_target'
                    if ~isnan(baseline_target_val)
                        baseline_val = baseline_target_val;
                        baseline_type_used = "target_specific";
                    else
                        baseline_val = baseline_global_val;
                        baseline_type_used = "global";
                    end

                otherwise
                    error('Unrecognized baseline_mode.');
            end

            if isnan(baseline_val)
                continue;
            end

            % --------------------------------------------------------
            % Collect trials for the target
            % --------------------------------------------------------
            firing_rate = [];

            for set_id = 1:n_sets
                idx_trials = find(([data(set_id).Data(array).Interp.Target_ID] == target_id) & ...
                                  ([data(set_id).Data(array).Interp.Excluded] == 0));

                for j = 1:length(idx_trials)
                    M_spikes = data(set_id).Data(array).Interp(idx_trials(j)).Trial(:, channel);
                    firing_rate = [firing_rate, M_spikes ./ bin_size];
                end
            end

            if isempty(firing_rate)
                continue;
            end

            % --------------------------------------------------------
            % Smoothing
            % --------------------------------------------------------
            data_smoothed = smoothdata(mean(firing_rate, 2), 'gaussian', w);
            t = (0:numel(data_smoothed)-1) * bin_size;

            % limit the window to the actual duration
            this_win_start = max(win_start, t(1));
            this_win_end   = min(win_end,   t(end));

            idx_win = t >= this_win_start & t <= this_win_end;
            if ~any(idx_win)
                continue;
            end

            t_win = t(idx_win);
            y_win = data_smoothed(idx_win);

            % --------------------------------------------------------
            % Response analysis using the function
            % --------------------------------------------------------
            response = find_response_peaks(t_win, y_win, baseline_val, min_delta_from_baseline, min_duration_bins);
            if ~(response.show_exc || response.show_inh)
                continue;
            end

            response_type = response.response_type;

            % Save in the target-specific struct
            responsive_channels(t_id).channels_global(end+1,1) = ch_global;
            responsive_channels(t_id).channels_local(end+1,1) = channel;
            responsive_channels(t_id).array_name(end+1,1) = array_names(array);
            responsive_channels(t_id).response_type(end+1,1) = response_type;
            responsive_channels(t_id).baseline_type(end+1,1) = baseline_type_used;

            % Save in a logic matrix
            responsive_matrix(t_id, ch_global) = true;

            % Save in a table
            new_row = table( ...
                target_id, ...
                string(array_names(array)), ...
                channel, ...
                ch_global, ...
                string(response_type), ...
                string(baseline_type_used), ...
                baseline_val, ...
                response.peak_exc_abs, response.peak_exc_rel, response.peak_exc_time, ...
                response.peak_inh_abs, response.peak_inh_rel, response.peak_inh_time, ...
                'VariableNames', { ...
                'target_id', 'array_name', 'channel_local', 'channel_global', ...
                'response_type', 'baseline_type', 'baseline', ...
                'peak_exc_abs', 'peak_exc_rel', 'peak_exc_time', ...
                'peak_inh_abs', 'peak_inh_rel', 'peak_inh_time'});

            results = [results; new_row]; 

            fprintf('Target %d | %s | ch %d (glob %d) | %s | baseline: %s\n', ...
                target_id, array_names(array), channel, ch_global, response_type, baseline_type_used);

        end
    end
end

%% Summary
fprintf('\n\n=========================================\n');
fprintf('SUMMARY OF RESPONSIVE CHANNELS BY TARGET\n');
fprintf('=========================================\n');

for t_id = 1:numel(target_des)
    fprintf('Target %d: %d responsive channels\n', ...
        responsive_channels(t_id).target_id, ...
        numel(responsive_channels(t_id).channels_global));
end

%% Save results
% To be changed based on the experimental condition:
save('responsive_channels_gaze_on_target_.mat', ...                                       
    'responsive_channels', 'results', 'responsive_matrix', ...
    'min_delta_from_baseline', 'min_duration_s', 'response_window_mode', ...
    'baseline_mode');