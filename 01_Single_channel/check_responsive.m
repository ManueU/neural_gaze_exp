clearvars -except mean_baseline_common mean_baseline_by_targets
close all
clc

% ============================================================
% FILE
% ============================================================
results_file = 'responsive_channels_gaze_on_center.mat';
data_file    = '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat';

load(results_file, 'responsive_matrix', 'results', ...
    'min_delta_from_baseline', 'min_duration_s', 'response_window_mode', 'baseline_mode');

S = load(data_file);
data = S.data;

% ============================================================
% PARAMETRI
% ============================================================
sets_plot = [1, 2, 3, 4, 5, 6];
n_sets = numel(sets_plot);

n_arrays   = 2;
n_channels = 96;
target_des = 1:8;
bin_size   = 0.02;
w          = 15;

array_names = ["medial", "lateral"];
plot_trial_mean = true;   % true = mostra anche PSTH medio non smussato
show_thresholds = true;

% Se vuoi limitare i plot a un sottoinsieme:
array_list   = 2:n_arrays;
channel_list = 83:n_channels;

% ============================================================
% EVENTI TASK
% ============================================================
TS = data(1).Data(1).Interp(1).Task_states;
state_names = string(TS(:,1));
state_dur_s = cellfun(@(x) size(x,1) * bin_size, TS(:,2));
state_onset_s = [0; cumsum(state_dur_s(1:end-1))];
get_onset = @(name) state_onset_s(find(strcmpi(state_names, name), 1, 'first'));

lines_t   = [get_onset("pres12"); get_onset("reach")];
lines_lab = ["Target cue"; "Go cue"];

idx_go = find(strcmp(lines_lab, "Go cue"), 1);
if isempty(idx_go)
    error('Go cue non trovato.');
end
go_time = lines_t(idx_go);

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
        error('response_window_mode non riconosciuto.');
end

% ============================================================
% LOOP
% ============================================================
for array = array_list
    for channel = channel_list

        ch_global = (array - 1) * n_channels + channel;

        for t_id = 1:numel(target_des)
            target_id = target_des(t_id);

            % -----------------------------------------------
            % Raccogli trial
            % -----------------------------------------------
            firing_rate = [];

            for set_ = 1:n_sets
                set_id = sets_plot(set_);
                idx_trials = find(([data(set_id).Data(array).Interp.Target_ID] == target_id) & ...
                                  ([data(set_id).Data(array).Interp.Excluded] == 0));

                for j = 1:length(idx_trials)
                    M_spikes = data(set_id).Data(array).Interp(idx_trials(j)).Trial(:, channel);
                    firing_rate = [firing_rate, M_spikes ./ bin_size]; 
                end
            end

            if isempty(firing_rate)
                fprintf('Array %s | ch %d (glob %d) | target %d -> nessun dato\n', ...
                    array_names(array), channel, ch_global, target_id);
                continue;
            end

            mean_fr = mean(firing_rate, 2);
            psth    = smoothdata(mean_fr, 'gaussian', w);
            t       = (0:numel(psth)-1) * bin_size;

            % -----------------------------------------------
            % Responsive o no?
            % -----------------------------------------------
            is_responsive = false;
            if target_id <= size(responsive_matrix,1) && ch_global <= size(responsive_matrix,2)
                is_responsive = logical(responsive_matrix(target_id, ch_global));
            end

            % Riga corrispondente in results
            row_idx = [];
            if exist('results', 'var') && ~isempty(results)
                row_idx = find(results.target_id == target_id & results.channel_global == ch_global, 1, 'first');
            end

            has_result_row = ~isempty(row_idx);

            % -----------------------------------------------
            % Baseline da mostrare
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
                % fallback su baseline globale se disponibile
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

            % --- subplot 1: PSTH completo
            subplot(1,2,1); hold on

            if plot_trial_mean
                plot(t, mean_fr, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
            end
            plot(t, psth, 'k', 'LineWidth', 1.8);

            % baseline
            if ~isnan(baseline_val)
                yline(baseline_val, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);
                if show_thresholds
                    yline(baseline_val + min_delta_from_baseline, ':', 'Color', 'r', 'LineWidth', 1.2);
                    yline(baseline_val - min_delta_from_baseline, ':', 'Color', 'r', 'LineWidth', 1.2);
                end
            end

            % eventi
            for k = 1:numel(lines_t)
                xline(lines_t(k), '--k', 'HandleVisibility', 'off');
            end

            % finestra analisi
            yl = ylim;
            patch([win_start win_end win_end win_start], [yl(1) yl(1) yl(2) yl(2)], ...
                  [0.85 0.9 1], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            if plot_trial_mean
                plot(t, mean_fr, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
            end
            plot(t, psth, 'k', 'LineWidth', 1.8);

            % picchi
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
            title('PSTH completo');
            grid on
            hold off

            % --- subplot 2: zoom finestra analisi
            subplot(1,2,2); hold on

            idx_win = t >= win_start & t <= win_end;
            if any(idx_win)
                if plot_trial_mean
                    plot(t(idx_win), mean_fr(idx_win), 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
                end
                plot(t(idx_win), psth(idx_win), 'k', 'LineWidth', 1.8);
            end

            if ~isnan(baseline_val)
                yline(baseline_val, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);
                if show_thresholds
                    yline(baseline_val + min_delta_from_baseline, ':', 'Color', 'r', 'LineWidth', 1.2);
                    yline(baseline_val - min_delta_from_baseline, ':', 'Color', 'r', 'LineWidth', 1.2);
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
            title('Zoom finestra analisi');
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
           

            % -----------------------------------------------
            % Console
            % -----------------------------------------------
            fprintf('\nArray %s | ch %d (glob %d) | target %d\n', ...
                array_names(array), channel, ch_global, target_id);
            fprintf('responsive_matrix = %d\n', is_responsive);

            if has_result_row
                fprintf('response_type = %s\n', response_type);
                fprintf('baseline = %.3f\n', baseline_val);
                fprintf('peak_exc_rel = %.3f | peak_exc_time = %.3f\n', peak_exc_rel, peak_exc_time);
                fprintf('peak_inh_rel = %.3f | peak_inh_time = %.3f\n', peak_inh_rel, peak_inh_time);
            else
                fprintf('nessuna riga corrispondente in results\n');
            end

            % -----------------------------------------------
            % Pausa manuale
            % -----------------------------------------------
            str = input(sprintf(['Plot completato.\n' ...
                                 'Premi INVIO per continuare, "q" per uscire, "s" per saltare al prossimo canale: ']), 's');

            if strcmpi(str, 'q')
                disp('Interruzione richiesta dall''utente.');
                return;
            elseif strcmpi(str, 's')
                break; % esce dal loop target e passa al prossimo canale
            end
        end
    end
end