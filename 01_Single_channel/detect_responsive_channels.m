clearvars -except mean_baseline_common std_baseline_common mean_baseline std_baseline
close all
clc

% ============================================================
% PARAMETRI
% ============================================================
sets_plot = [1, 2, 3, 4, 5, 6];
n_sets = numel(sets_plot);

n_arrays = 2;
n_channels = 96;
ch_start = 1;
ch_end = 96;

target_des = 1:8;
bin_size = 0.02;

filename = '../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat'; %% DA MODIFICARE IN BASE ALLA CONDIZIONE

array_names = ["medial", "lateral"];

% smoothing
w = 15;

% finestra di ricerca nella condizione free-gaze
% 'planning'  = Go cue - 0.5 s --> Go cue
% 'reaching'  = Go cue         --> Go cue + 3.0 s
% 'full'      = Go cue - 0.5 s --> Go cue + 3.0 s
response_window_mode = 'full';

% soglia assoluta rispetto alla baseline
min_delta_from_baseline = 12;   % Hz

% durata minima della modulazione
min_duration_s = 0.06;          % 60 ms
min_duration_bins = ceil(min_duration_s / bin_size);

% baseline mode:
% 'global'          -> usa mean_baseline_common
% 'target_specific' -> usa mean_baseline.by_targets(target, ch_global)
% 'prefer_target'   -> usa target-specifica se disponibile, altrimenti globale
baseline_mode = 'global';


% ============================================================
% CHECK BASELINE
% ============================================================
has_global_baseline = exist('mean_baseline_common', 'var') && ~isempty(mean_baseline_common);
has_target_baseline = exist('mean_baseline', 'var') && ~isempty(mean_baseline.by_targets);

if ~has_global_baseline && ~has_target_baseline
    error('Nessuna baseline disponibile: serve mean_baseline_common e/o mean_baseline_by_targets.');
end

% ============================================================
% CARICAMENTO DATI
% ============================================================
S = load(filename);
data = S.data;

TS = data(1).Data(1).Interp(1).Task_states;
state_names = string(TS(:,1));
state_dur_s = cellfun(@(x) size(x,1) * bin_size, TS(:,2));

state_onset_s = [0; cumsum(state_dur_s(1:end-1))];
get_onset = @(name) state_onset_s(find(strcmpi(state_names, name), 1, 'first'));


task_name = "Free-gaze";                                       %% DA MODIFICARE IN BASE ALLA CONDIZIONE
lines_t   = [get_onset("pres12"); get_onset("reach")];
lines_lab = ["Target cue"; "Go cue"];

idx_go = find(strcmp(lines_lab, "Go cue"), 1);
if isempty(idx_go)
    error('Go cue non trovato nei Task_states del file gaze-on-target.');
end
go_time = lines_t(idx_go);

% definizione finestra
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
% STRUTTURE OUTPUT
% ============================================================
responsive_channels = struct();

results = table();
responsive_matrix = false(numel(target_des), n_arrays*n_channels);

% ============================================================
% LOOP TARGET / ARRAY / CANALE
% ============================================================
for t_id = 1:numel(target_des)
    target_id = target_des(t_id);

    responsive_channels(t_id).target_id = target_id;
    responsive_channels(t_id).channels_global = [];
    responsive_channels(t_id).channels_local = [];
    responsive_channels(t_id).array_name = strings(0,1);
    responsive_channels(t_id).response_type = strings(0,1);
    responsive_channels(t_id).baseline_type = strings(0,1);

    fprintf('\n==============================\n');
    fprintf('TARGET %d\n', target_id);
    fprintf('==============================\n');

    for array = 1:n_arrays
        for channel = ch_start:ch_end

            ch_global = (array - 1) * n_channels + channel;

            % --------------------------------------------------------
            % Selezione baseline
            % --------------------------------------------------------
            baseline_val = NaN;
            baseline_type_used = "";

            baseline_global_val = NaN;
            if has_global_baseline && ch_global <= numel(mean_baseline_common)
                baseline_global_val = mean_baseline_common(ch_global);
            end

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
                    error('baseline_mode non riconosciuto.');
            end

            if isnan(baseline_val)
                continue;
            end

            % --------------------------------------------------------
            % Raccogli trial del target
            % --------------------------------------------------------
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
                continue;
            end

            % --------------------------------------------------------
            % PSTH medio smussato
            % --------------------------------------------------------
            data_smoothed = smoothdata(mean(firing_rate, 2), 'gaussian', w);
            t = (0:numel(data_smoothed)-1) * bin_size;

            % limita la finestra alla durata reale
            this_win_start = max(win_start, t(1));
            this_win_end   = min(win_end,   t(end));

            idx_win = t >= this_win_start & t <= this_win_end;
            if ~any(idx_win)
                continue;
            end

            t_win = t(idx_win);
            y_win = data_smoothed(idx_win);

            % --------------------------------------------------------
            % Analisi risposta tramite funzione
            % --------------------------------------------------------
            response = find_response_peaks( ...
                t_win, y_win, baseline_val, ...
                min_delta_from_baseline, min_duration_bins);

            if ~(response.show_exc || response.show_inh)
                continue;
            end

            % tipo risposta
            response_type = response.response_type;

            % salva in struct per target
            responsive_channels(t_id).channels_global(end+1,1) = ch_global;
            responsive_channels(t_id).channels_local(end+1,1) = channel;
            responsive_channels(t_id).array_name(end+1,1) = array_names(array);
            responsive_channels(t_id).response_type(end+1,1) = response_type;
            responsive_channels(t_id).baseline_type(end+1,1) = baseline_type_used;

            % salva in matrice logica
            responsive_matrix(t_id, ch_global) = true;

            % salva in tabella lunga
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

% ============================================================
% RIEPILOGO
% ============================================================
fprintf('\n\n=========================================\n');
fprintf('RIEPILOGO NUMERO CANALI RESPONSIVI PER TARGET\n');
fprintf('=========================================\n');

for t_id = 1:numel(target_des)
    fprintf('Target %d: %d canali responsivi\n', ...
        responsive_channels(t_id).target_id, ...
        numel(responsive_channels(t_id).channels_global));
end

% ============================================================
% SALVATAGGIO
% ============================================================
save('responsive_channels_free_gaze_prova.mat', ...                                       %% DA MODIFICARE IN BASE ALLA CONDIZIONE
    'responsive_channels', 'results', 'responsive_matrix', ...
    'min_delta_from_baseline', 'min_duration_s', 'response_window_mode', ...
    'baseline_mode');

disp('File salvato: responsive_channels_free_gaze_prova.mat');

% ============================================================
% HEATMAP
% ============================================================
figure('Color','w');
imagesc(responsive_matrix);
colormap(gray);
xlabel('Channel');
ylabel('Target');
title(sprintf('Responsive channels in free-gaze'));
yticks(1:numel(target_des));
yticklabels(string(target_des));