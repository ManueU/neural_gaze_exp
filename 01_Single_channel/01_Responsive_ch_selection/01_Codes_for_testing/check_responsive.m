%% Responsive channel inspection

% Interactive visualization of PSTHs, response peaks, and baseline thresholds
% for channels classified as responsive during the detection analysis.

clearvars -except mean_baseline_common mean_baseline_by_targets
close all
clc

%% File selection and data loading
% To be changed based on the experimental condition:
results_file = 'responsive_channels_free-gaze.mat';
% results_file = 'responsive_channels_gaze_on_center.mat';
% results_file = 'responsive_channels_gaze_on_target.mat';

data_file    = '../../../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat';
% data_file    = '../../../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat';
% data_file    = '../../../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat';

load(results_file, 'responsive_matrix', 'results', ...
    'min_delta_from_baseline', 'min_duration_s', 'response_window_mode', 'baseline_mode');

S = load(data_file);
data = S.data;

%% Parameters
n_sets = 6;
n_arrays   = 2;
n_channels = 96;
target_des = 1:8;
bin_size   = 0.02;
w          = 15;

array_names = ["medial", "lateral"];
plot_trial_mean = true;   % true = also display the unsmoothed mean PSTH
show_thresholds = true;

array_list   = 2:n_arrays;
channel_list = 1:n_channels;

%% Task events
TS = data(1).Data(1).Interp(1).Task_states;
state_names = string(TS(:,1));
state_dur_s = cellfun(@(x) size(x,1) * bin_size, TS(:,2));
state_onset_s = [0; cumsum(state_dur_s(1:end-1))];
get_onset = @(name) state_onset_s(find(strcmpi(state_names, name), 1, 'first'));

% To be changed based on the experimental condition:
task_name = "Free-gaze"; 
lines_t   = [get_onset("pres12"); get_onset("reach")];
lines_lab = ["Target cue"; "Go cue"];
%------
% task_name = "Gaze-on-center"; 
% lines_t   = [get_onset("pres12"); get_onset("reach")];
% lines_lab = ["Target cue"; "Go cue"];
%------
% task_name = "Gaze-on-target"; 
% lines_t   = [get_onset("Gaze"); get_onset("reach")];
% lines_lab = ["Target cue"; "Go cue"];

idx_go = find(strcmp(lines_lab, "Go cue"), 1);
if isempty(idx_go)
    error('Go cue non trovato.');
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

%% Target / Array / Channel loop
for array = array_list
    for channel = channel_list

        ch_global = (array - 1) * n_channels + channel;

        for t_id = 1:numel(target_des)
            target_id = target_des(t_id);

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
                fprintf('Array %s | ch %d (glob %d) | target %d -> no data\n', ...
                    array_names(array), channel, ch_global, target_id);
                continue;
            end

            mean_fr = mean(firing_rate, 2);
            psth    = smoothdata(mean_fr, 'gaussian', w);
            t       = (0:numel(psth)-1) * bin_size;

            % -----------------------------------------------
            % Responsive or not?
            % -----------------------------------------------
            is_responsive = false;
            if target_id <= size(responsive_matrix,1) && ch_global <= size(responsive_matrix,2)
                is_responsive = logical(responsive_matrix(target_id, ch_global));
            end

            row_idx = [];
            if exist('results', 'var') && ~isempty(results)
                row_idx = find(results.target_id == target_id & results.channel_global == ch_global, 1, 'first');
            end

            has_result_row = ~isempty(row_idx);

            % -----------------------------------------------
            % Baseline
            % -----------------------------------------------
            baseline_val = NaN;
            response_type = "none";
            peak_exc_abs = NaN; peak_exc_rel = NaN; peak_exc_time = NaN;
            peak_inh_abs = NaN; peak_inh_rel = NaN; peak_inh_time = NaN;
            baseline_type_used = "";

            if has_result_row
                baseline_val       = results.baseline(row_idx);
                response_type      = string(results.response_type(row_idx));
                peak_exc_abs       = results.peak_exc_abs(row_idx);
                peak_exc_rel       = results.peak_exc_rel(row_idx);
                peak_exc_time      = results.peak_exc_time(row_idx);
                peak_inh_abs       = results.peak_inh_abs(row_idx);
                peak_inh_rel       = results.peak_inh_rel(row_idx);
                peak_inh_time      = results.peak_inh_time(row_idx);

                if any(strcmp('baseline_type', results.Properties.VariableNames))
                    baseline_type_used = string(results.baseline_type(row_idx));
                end
            else
                if exist('mean_baseline_common', 'var') && ~isempty(mean_baseline_common) && ch_global <= numel(mean_baseline_common)
                    baseline_val = mean_baseline_common(ch_global);
                    baseline_type_used = "global-fallback";
                end
            end

            % -----------------------------------------------
            % Plot
            % -----------------------------------------------
            close all
            figure('Color', 'w', 'Position', [100 100 1200 500]);

            % --- subplot 1: PSTH
            subplot(1,2,1); hold on

            if plot_trial_mean
                plot(t, mean_fr, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
            end
            plot(t, psth, 'k', 'LineWidth', 1.8);

            % baseline
            if ~isnan(baseline_val)
                yline(baseline_val, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
                if show_thresholds
                    yline(baseline_val + min_delta_from_baseline, ':', 'Color', 'r', 'LineWidth', 2);
                    yline(baseline_val - min_delta_from_baseline, ':', 'Color', 'r', 'LineWidth', 2);
                end
            end

            % events
            for k = 1:numel(lines_t)
                xline(lines_t(k), '--k', 'HandleVisibility', 'off');
            end

            % window of analysis
            yl = ylim;
            patch([win_start win_end win_end win_start], [yl(1) yl(1) yl(2) yl(2)], ...
                  [0.85 0.9 1], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            if plot_trial_mean
                plot(t, mean_fr, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
            end
            plot(t, psth, 'k', 'LineWidth', 1.8);

            % peaks
            if has_result_row
                if ~isnan(peak_exc_abs)
                    scatter(peak_exc_time, peak_exc_abs, 70, 'r', 'filled');
                end
                if ~isnan(peak_inh_abs)
                    scatter(peak_inh_time, peak_inh_abs, 70, 'b', 'filled');
                end
            end

            xlabel('Time (s)');
            ylabel('Firing rate (Hz)');
            title('Full PSTH');
            grid on
            hold off

            % --- subplot 2: zoom
            subplot(1,2,2); hold on

            idx_win = t >= win_start & t <= win_end;
            if any(idx_win)
                if plot_trial_mean
                    plot(t(idx_win), mean_fr(idx_win), 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
                end
                plot(t(idx_win), psth(idx_win), 'k', 'LineWidth', 1.8);
            end

            if ~isnan(baseline_val)
                yline(baseline_val, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
                if show_thresholds
                    yline(baseline_val + min_delta_from_baseline, ':', 'Color', 'r', 'LineWidth', 2);
                    yline(baseline_val - min_delta_from_baseline, ':', 'Color', 'r', 'LineWidth', 2);
                end
            end

            if has_result_row
                if ~isnan(peak_exc_abs)
                    scatter(peak_exc_time, peak_exc_abs, 70, 'r', 'filled');
                end
                if ~isnan(peak_inh_abs)
                    scatter(peak_inh_time, peak_inh_abs, 70, 'b', 'filled');
                end
            end

            xlabel('Time (s)');
            ylabel('Firing rate (Hz)');
            title('Zoom');
            grid on
            hold off

            % -----------------------------------------------
            % Titolo generale
            % -----------------------------------------------
            if is_responsive
                status_txt = 'Responsive';
            else
                status_txt = 'Non responsive';
            end

            sgtitle(sprintf('%s | ch %d (glob %d) | target %d | %s', ...
                                 array_names(array), channel, ch_global, target_id, status_txt), ...
                                 'FontWeight', 'bold');
          
            str = input(sprintf(['Plot completed.\n' ...
                    'Press ENTER to continue, "q" to quit, or "s" to skip to the next channel: ']), 's');

            if strcmpi(str, 'q')
                disp('Execution interrupted by the user');
                return;
            elseif strcmpi(str, 's')
                break; % Exit the target loop and move to the next channel
            end
        end
    end
end