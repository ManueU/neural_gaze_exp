clearvars -except mean_baseline_common std_baseline_common mean_baseline std_baseline
close all
clc

sets_plot = [1, 2, 3, 4, 5, 6];
n_sets = numel(sets_plot);
n_arrays = 2;
n_channels = 96;
ch_start = 1;
ch_end = 96;
target_des = [1,2,3,4,5,6,7,8];
bin_size = 0.02;

filename = { ...
    '../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat', ...
    '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat', ...
    '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat'};

array_names = ["medial", "lateral"];
colors_target = [
    0.839, 0.153, 0.157;  % rosso
    0.122, 0.467, 0.706;  % blu
    0.172, 0.627, 0.172;  % verde
    0.580, 0.404, 0.741;  % viola
    1.000, 0.498, 0.055;  % arancione
    0.737, 0.741, 0.133;  % giallo oliva
    0.549, 0.337, 0.294;  % marrone
    0.890, 0.466, 0.760;  % rosa
];

% ------------------------------------------------------------
% BASELINE MODE
% 'global'          -> usa mean_baseline_common
% 'target_specific' -> usa mean_baseline_by_targets(target, ch_global)
% 'prefer_target'   -> usa target-specifica se disponibile, altrimenti globale
% ------------------------------------------------------------
baseline_mode = 'global';

% Se true, gli zeri nella matrice target-specifica vengono trattati come mancanti
treat_zero_target_baseline_as_missing = false;

has_global_baseline = exist('mean_baseline_common', 'var') && ~isempty(mean_baseline_common);
has_target_baseline = exist('mean_baseline', 'var') && ~isempty(mean_baseline.by_targets);

if ~has_global_baseline && ~has_target_baseline
    error('Nessuna baseline disponibile: serve mean_baseline_common e/o mean_baseline_by_targets.');
end

% ------------------------------------------------------------
% OPZIONE: scegliere tutti i canali oppure un solo canale globale
% ------------------------------------------------------------
use_single_global_channel = true;   % true = mostra solo un canale globale
selected_global_channel   = 169;     % canale globale da mostrare se true

if use_single_global_channel
    if selected_global_channel < 1 || selected_global_channel > n_arrays*n_channels
        error('selected_global_channel deve essere compreso tra 1 e %d', n_arrays*n_channels);
    end

    array_list = ceil(selected_global_channel / n_channels);
    channel_list_per_array = selected_global_channel - (array_list - 1) * n_channels;
else
    array_list = 1:n_arrays;
    channel_list_per_array = ch_start:ch_end;
end

% ------------------------------------------------------------
% Caricamento di tutte le condizioni
% ------------------------------------------------------------
cond_data = struct();

for d = 1:numel(filename)
    S = load(filename{d});
    data = S.data;

    TS = data(1).Data(1).Interp(1).Task_states;
    state_names = string(TS(:,1));
    state_dur_s = cellfun(@(x) size(x,1) * bin_size, TS(:,2));

    state_onset_s = [0; cumsum(state_dur_s(1:end-1))];
    get_onset = @(name) state_onset_s(find(strcmpi(state_names, name), 1, 'first'));

    thisfile = string(filename{d});

    if contains(thisfile, "free-gaze", "IgnoreCase", true)
        task_name = "Free-gaze";
        lines_t   = [get_onset("pres12"); get_onset("reach")];
        lines_lab = ["Target cue"; "Go cue"];

    elseif contains(thisfile, "motor", "IgnoreCase", true)
        task_name = "Gaze-on-center";
        lines_t   = [get_onset("pres12"); get_onset("reach")];
        lines_lab = ["Target cue"; "Go cue"];

    elseif contains(thisfile, "controlled", "IgnoreCase", true)
        task_name = "Gaze-on-target";
        lines_t   = [get_onset("pres12"); get_onset("gaze"); get_onset("reach")];
        lines_lab = ["Target cue"; "Go cue - gaze"; "Go cue"];

    elseif contains(thisfile, "gaze", "IgnoreCase", true)
        task_name = "Gaze-only";
        lines_t   = [get_onset("pres12"); get_onset("gaze")];
        lines_lab = ["Target cue"; "Go cue - gaze"];
    else
        task_name = "Unknown";
        lines_t   = [];
        lines_lab = strings(0,1);
    end

    cond_data(d).data = data;
    cond_data(d).task_name = task_name;
    cond_data(d).lines_t = lines_t;
    cond_data(d).lines_lab = lines_lab;
end

% ------------------------------------------------------------
% Loop principale
% Una figure per ogni combinazione canale-target
% Tre subplot = tre condizioni
% ------------------------------------------------------------
w = 15;

for array = array_list

    if use_single_global_channel
        channel_values = channel_list_per_array;
    else
        channel_values = ch_start:ch_end;
    end

    for channel = channel_values

        ch_global = (array - 1) * n_channels + channel;

        fprintf('\nArray: %s | Channel locale: %d | Channel globale: %d\n', ...
            array_names(array), channel, ch_global);

        for target = 1:length(target_des)

            target_id = target_des(target);

            close all;
            figure('Color', 'w', ...
                   'Name', sprintf('Array %s | Ch %d | Target %d', array_names(array), ch_global, target_id), ...
                   'Position', [100 100 1400 450]);

            sgtitle(sprintf('Array: %s | Channel %d | Target %d', ...
                array_names(array), ch_global, target_id), ...
                'FontWeight', 'bold');

            y_limits_all = nan(numel(cond_data), 2);

            for d = 1:numel(cond_data)

                % --------------------------------------------------------
                % Baseline del canale corrente
                % --------------------------------------------------------
                baseline_val = NaN;
                baseline_type_used = "";
    
                baseline_global_val = NaN;
                if has_global_baseline && ch_global <= numel(mean_baseline_common)
                    baseline_global_val = mean_baseline_common(ch_global);
                end
    
                baseline_target_val = NaN;
                if has_target_baseline
                    if target_id <= size(mean_baseline.by_targets{1,d},1) && ch_global <= size(mean_baseline.by_targets{1,d},2)
                        baseline_target_val = mean_baseline.by_targets{1,d}(target_id, ch_global);
    
                        if treat_zero_target_baseline_as_missing && baseline_target_val == 0
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
                        baseline_type_used = "target specific";
    
                    case 'prefer_target'
                        if ~isnan(baseline_target_val)
                            baseline_val = baseline_target_val;
                            baseline_type_used = "target specific";
                        else
                            baseline_val = baseline_global_val;
                            baseline_type_used = "global";
                        end
    
                    otherwise
                        error('baseline_mode non riconosciuto.');
                end

                data = cond_data(d).data;
                task_name = cond_data(d).task_name;
                lines_t = cond_data(d).lines_t;
                lines_lab = cond_data(d).lines_lab;

                firing_rate = [];

                for set_ = 1:n_sets
                    set_id = sets_plot(set_);
                    idx = find(([data(set_id).Data(array).Interp.Target_ID] == target_id) & ...
                               ([data(set_id).Data(array).Interp.Excluded] == 0));

                    for j = 1:length(idx)
                        M_spikes = data(set_id).Data(array).Interp(idx(j)).Trial(:, channel);
                        firing_rate = [firing_rate, M_spikes ./ bin_size]; 
                    end
                end

                subplot(1, numel(cond_data), d);
                hold on;

                if ~isempty(firing_rate)
                    data_smoothed = smoothdata(mean(firing_rate, 2), 'gaussian', w);
                    t = (0:numel(data_smoothed)-1) * bin_size;

                    plot(t, data_smoothed, ...
                        'LineWidth', 1.8, ...
                        'Color', colors_target(target, :));

                    % Linea orizzontale di baseline
                    yline(baseline_val, '--', ...
                        sprintf('Baseline'), ...
                        'Color', [0.2 0.2 0.2], ...
                        'LineWidth', 1.2, ...
                        'HandleVisibility', 'off');
                else
                    text(0.5, 0.5, 'Nessun dato', ...
                        'HorizontalAlignment', 'center', ...
                        'Units', 'normalized', ...
                        'FontSize', 12);
                    t = [];
                end

                for k = 1:numel(lines_t)
                    xline(lines_t(k), '--k', 'HandleVisibility', 'off');

                    if k <= numel(lines_lab) && strcmp(lines_lab(k), "Go cue - gaze")
                        xline(lines_t(k) + 0.08, '--k', 'HandleVisibility', 'off');
                    end
                     
                    % if k <= numel(lines_lab) && strcmp(lines_lab(k), "Go cue")
                    %     xline(lines_t(k) - 0.5, '-k', 'HandleVisibility', 'off');
                    %     xline(lines_t(k) + 0.5, '-k', 'HandleVisibility', 'off');
                    % end
                end

                title(char(task_name), 'Interpreter', 'none');
                xlabel('Time (s)');
                ylabel('Firing rate');
                axis tight;

                ax = gca;

                if ~isempty(t)
                    y_limits_all(d, :) = ax.YLim;
                end

                hold off;
            end

            % Uniforma l'asse Y nei 3 subplot per confronto diretto
            valid_rows = all(~isnan(y_limits_all), 2);
            if any(valid_rows)
                y_min = min(y_limits_all(valid_rows, 1));
                y_max = max(y_limits_all(valid_rows, 2));

                if y_min == y_max
                    y_min = y_min - 1;
                    y_max = y_max + 1;
                end

                for d = 1:numel(cond_data)
                    subplot(1, numel(cond_data), d);
                    ylim([y_min y_max]);
                end

                % Banda grigia tra -500 ms e +500 ms dal Go cue
                % for d = 1:numel(cond_data)
                %     subplot(1, numel(cond_data), d);
                % 
                %     lines_t = cond_data(d).lines_t;
                %     lines_lab = cond_data(d).lines_lab;
                % 
                %     idx_go = find(strcmp(lines_lab, "Go cue"), 1);
                % 
                %     if ~isempty(idx_go)
                %         xl1 = lines_t(idx_go) - 0.5;
                %         xl2 = lines_t(idx_go) + 0.5;
                % 
                %         p = patch([xl1 xl2 xl2 xl1], [y_min y_min y_max y_max], [0.8 0.8 0.8], ...
                %                   'FaceAlpha', 0.3, ...
                %                   'EdgeColor', 'none', ...
                %                   'HandleVisibility', 'off');
                % 
                %         uistack(p, 'bottom');
                %     end
                % end

                % ------------------------------------------------------------
                % Etichette delle fasi: stessa coordinata Y per tutti i subplot
                % ------------------------------------------------------------
                y_txt = y_max - 0.05 * (y_max - y_min);

                for d = 1:numel(cond_data)
                    subplot(1, numel(cond_data), d);
                    ax = gca;

                    lines_t = cond_data(d).lines_t;
                    lines_lab = cond_data(d).lines_lab;

                    if ~isempty(lines_t)
                        dx = 0.06 * max(range(ax.XLim), eps);

                        text(lines_t - dx, y_txt * ones(size(lines_t)), lines_lab, ...
                            'HorizontalAlignment', 'right', ...
                            'VerticalAlignment', 'top', ...
                            'Rotation', 90, ...
                            'Clipping', 'on', ...
                            'HandleVisibility', 'off');
                    end
                end
            end

            str = input(sprintf(['Canale %d - Target %d completato.\n' ...
                                 'Baseline usata: %s\n' ...
                                 'Premi INVIO per continuare, oppure digita q e premi INVIO per uscire: '], ...
                                 ch_global, target_id, baseline_type_used), 's');

            if strcmpi(str, 'q')
                disp('Interruzione richiesta dall''utente.');
                return;
            end
        end
    end
end