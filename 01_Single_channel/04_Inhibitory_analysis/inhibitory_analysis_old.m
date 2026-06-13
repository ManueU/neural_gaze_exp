% clearvars -except mean_baseline_common
% close all
clc

% ============================================================
% LOAD
% ============================================================
load('responsive_channels_free_gaze.mat', 'responsive_channels', 'responsive_matrix');

if ~exist('responsive_channels','var') || isempty(responsive_channels)
    error('responsive_channels non disponibile.');
end

if ~exist('mean_baseline_common','var') || isempty(mean_baseline_common)
    error('Serve mean_baseline_common nel workspace.');
end

% ============================================================
% FILE
% ============================================================
filename = '../00_Data_extraction/BCI02_Session_0924/gaze_BCI02_exclUpdated.mat';

% ============================================================
% PARAMETRI
% ============================================================
sets_plot = [1, 2, 3, 4, 5, 6];
n_sets = numel(sets_plot);

n_arrays = 2;
n_channels = 96;
target_des = 1:8;
bin_size = 0.02;
w = 15;

array_names = ["medial", "lateral"];

% finestra rispetto all'onset di "gaze"
win_start_rel = -1.0;   % s
win_end_rel   =  1.0;   % s

% criteri di inibizione
min_delta_from_baseline = 12;   % Hz
min_duration_s = 0.06;          % s
min_duration_bins = ceil(min_duration_s / bin_size);

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

gaze_time = get_onset("gaze");

if isempty(gaze_time) || isnan(gaze_time)
    error('Fase "gaze" non trovata nei Task_states.');
end

fprintf('\n============================================\n');
fprintf('GAZE-ONLY INHIBITORY ANALYSIS\n');
fprintf('============================================\n');
fprintf('Finestra: [%.2f %.2f] s rispetto a gaze onset\n', win_start_rel, win_end_rel);
fprintf('Soglia inibitoria: %.1f Hz sotto baseline\n', min_delta_from_baseline);
fprintf('Durata minima: %.0f ms (%d bins)\n', 1000*min_duration_s, min_duration_bins);

% ============================================================
% OUTPUT
% ============================================================
all_rows = table();

% ============================================================
% LOOP TARGET
% ============================================================
for t_id = 1:numel(target_des)

    target_id = target_des(t_id);

    sel_channels = responsive_channels(t_id).channels_global;

    if isempty(sel_channels)
        fprintf('\nTarget %d: nessun canale responsivo in free-gaze\n', target_id);
        continue;
    end

    fprintf('\n========================================\n');
    fprintf('TARGET %d | N canali responsivi free-gaze = %d\n', target_id, numel(sel_channels));
    fprintf('========================================\n');

    for c = 1:numel(sel_channels)

        ch_global = sel_channels(c);
        array = ceil(ch_global / n_channels);
        channel = ch_global - (array - 1) * n_channels;

        if ch_global > numel(mean_baseline_common)
            continue;
        end

        baseline_val = mean_baseline_common(ch_global);

        if isnan(baseline_val)
            continue;
        end

        % ----------------------------------------------------
        % Raccogli trial del target nella condizione gaze-only
        % ----------------------------------------------------
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

        % ----------------------------------------------------
        % PSTH medio smussato
        % ----------------------------------------------------
        psth = smoothdata(mean(firing_rate, 2), 'gaussian', w);
        t = (0:numel(psth)-1) * bin_size;

        this_win_start = gaze_time + win_start_rel;
        this_win_end   = gaze_time + win_end_rel;

        this_win_start = max(this_win_start, t(1));
        this_win_end   = min(this_win_end,   t(end));

        idx_win = t >= this_win_start & t <= this_win_end;
        if ~any(idx_win)
            continue;
        end

        t_win_abs = t(idx_win);
        y_win = psth(idx_win);
        t_win_rel = t_win_abs - gaze_time;

        % ----------------------------------------------------
        % Analisi inibitoria
        % ----------------------------------------------------
        response = find_inhibitory_response( ...
            t_win_rel, y_win, baseline_val, ...
            min_delta_from_baseline, min_duration_bins);

        new_row = table( ...
            target_id, ...
            string(array_names(array)), ...
            channel, ...
            ch_global, ...
            baseline_val, ...
            response.inh_present, ...
            response.peak_inh_abs, ...
            response.peak_inh_rel, ...
            response.peak_inh_amp, ...
            response.peak_inh_time, ...
            response.inh_area_total, ...
            response.inh_area_thresholded, ...
            response.n_valid_segments, ...
            'VariableNames', { ...
            'target_id', ...
            'array_name', ...
            'channel_local', ...
            'channel_global', ...
            'baseline', ...
            'inh_present', ...
            'peak_inh_abs', ...
            'peak_inh_rel', ...
            'peak_inh_amp', ...
            'peak_inh_time', ...
            'inh_area_total', ...
            'inh_area_thresholded', ...
            'n_valid_segments'});

        all_rows = [all_rows; new_row];

        if response.inh_present
            fprintf(['Target %d | %s | ch %d (glob %d) | INH | ' ...
                     'amp = %.2f Hz | time = %.2f s | area = %.2f Hz*s\n'], ...
                     target_id, array_names(array), channel, ch_global, ...
                     response.peak_inh_amp, response.peak_inh_time, response.inh_area_total);
        end
    end
end

% ============================================================
% SUMMARY PER TARGET
% ============================================================
fprintf('\n============================================\n');
fprintf('SUMMARY PER TARGET\n');
fprintf('============================================\n');

for t_id = 1:numel(target_des)
    Tt = all_rows(all_rows.target_id == target_des(t_id), :);

    if isempty(Tt)
        fprintf('Target %d: nessun dato\n', target_des(t_id));
        continue;
    end

    n_inh = sum(Tt.inh_present);
    n_tot = height(Tt);

    fprintf('Target %d | canali-target analizzati = %d | inibitori = %d (%.1f%%)\n', ...
        target_des(t_id), n_tot, n_inh, 100*n_inh/n_tot);
end

% ============================================================
% SUMMARY GLOBALE
% ============================================================
fprintf('\n============================================\n');
fprintf('SUMMARY GLOBALE\n');
fprintf('============================================\n');

n_total_rows = height(all_rows);
n_inh_rows = sum(all_rows.inh_present);

fprintf('Canali-target analizzati: %d\n', n_total_rows);
fprintf('Canali-target con inibizione: %d (%.1f%%)\n', ...
    n_inh_rows, 100*n_inh_rows/n_total_rows);

% canali regardless target
if ~isempty(all_rows)
    [G, channel_ids] = findgroups(all_rows.channel_global);
    ch_inh_any = splitapply(@(x) any(x), all_rows.inh_present, G);

    n_channels_inh_any = sum(ch_inh_any);
    fprintf('Canali con inibizione in almeno un target: %d\n', n_channels_inh_any);
end

% ============================================================
% SALVATAGGIO
% ============================================================
save('gaze_only_inhibitory_analysis_target_by_target.mat', ...
    'all_rows', ...
    'win_start_rel', 'win_end_rel', ...
    'min_delta_from_baseline', 'min_duration_s', ...
    'gaze_time');

disp('File salvato: gaze_only_inhibitory_analysis_target_by_target.mat');

% usa solo i casi con inibizione presente per ampiezza/area
T_inh = all_rows(all_rows.inh_present == true, :);

% ============================================================
% SUMMARY PER TARGET
% ============================================================
n_targets = numel(target_des);

pct_inh = nan(1, n_targets);
n_inh   = zeros(1, n_targets);
n_tot   = zeros(1, n_targets);

median_amp  = nan(1, n_targets);
median_area = nan(1, n_targets);

for i = 1:n_targets
    target_id = target_des(i);

    Tt = all_rows(all_rows.target_id == target_id, :);
    Tt_inh = T_inh(T_inh.target_id == target_id, :);

    n_tot(i) = height(Tt);
    n_inh(i) = height(Tt_inh);

    if n_tot(i) > 0
        pct_inh(i) = 100 * n_inh(i) / n_tot(i);
    end

    if ~isempty(Tt_inh)
        median_amp(i)  = median(Tt_inh.peak_inh_amp, 'omitnan');
        median_area(i) = median(Tt_inh.inh_area_total, 'omitnan');
    end
end

% ============================================================
% FIGURA
% ============================================================
figure('Color','w','Position',[100 100 1400 420]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% ============================================================
% 1) Percentuale di canali-target con inibizione per target
% ============================================================
nexttile
bar(target_des, pct_inh, 0.7);
xlabel('Target');
ylabel('Percentage of channels (%)');
title('Responsive channels by target');
xticks(target_des);
ylim([0 100]);
grid on
box off

for i = 1:n_targets
    if ~isnan(pct_inh(i))
        text(target_des(i), pct_inh(i)+3, sprintf('%.1f%%\n(n=%d/%d)', ...
            pct_inh(i), n_inh(i), n_tot(i)), ...
            'HorizontalAlignment','center', 'FontSize',9);
    end
end

% ============================================================
% 2) Boxplot ampiezza del picco inibitorio per target
% ============================================================
nexttile
if ~isempty(T_inh)
    boxplot(T_inh.peak_inh_amp, T_inh.target_id, 'Labels', string(target_des));
    hold on

    % mediana sovrapposta
    for i = 1:n_targets
        if ~isnan(median_amp(i))
            plot(i, median_amp(i), 'kd', 'MarkerFaceColor','k', 'MarkerSize',5);
        end
    end
    hold off
    
    ylim([0 80])
    xlabel('Target');
    ylabel('Peak amplitude (Hz)');
    title('Peak inhibitory amplitude');
    grid on
    box off
else
    text(0.5, 0.5, 'No inhibitory responses found', ...
        'Units','normalized', 'HorizontalAlignment','center');
    axis off
end

% ============================================================
% 3) Boxplot area sotto baseline per target
% ============================================================
nexttile
if ~isempty(T_inh)
    boxplot(T_inh.inh_area_total, T_inh.target_id, 'Labels', string(target_des));
    hold on

    % mediana sovrapposta
    for i = 1:n_targets
        if ~isnan(median_area(i))
            plot(i, median_area(i), 'kd', 'MarkerFaceColor','k', 'MarkerSize',5);
        end
    end
    hold off

    ylim([0 80])
    xlabel('Target');
    ylabel('Inhibitory area (Hz*s)');
    title('Area below baseline');
    grid on
    box off
else
    text(0.5, 0.5, 'No inhibitory responses found', ...
        'Units','normalized', 'HorizontalAlignment','center');
    axis off
end

sgtitle('Gaze-only characterization');


% ============================================================
% FUNZIONI LOCALI
% ============================================================
function response = find_inhibitory_response(t_win, y_win, baseline_val, min_delta, min_duration_bins)

    response = struct();
    response.inh_present = false;
    response.peak_inh_abs = NaN;
    response.peak_inh_rel = NaN;
    response.peak_inh_amp = NaN;
    response.peak_inh_time = NaN;
    response.inh_area_total = NaN;
    response.inh_area_thresholded = NaN;
    response.n_valid_segments = 0;

    % --------------------------------------------------------
    % Area sotto baseline nell'intera finestra
    % (solo porzioni sotto baseline)
    % --------------------------------------------------------
    below_baseline = max(baseline_val - y_win, 0);
    response.inh_area_total = trapz(t_win, below_baseline);

    % --------------------------------------------------------
    % Criterio di soglia: almeno min_delta sotto baseline
    % --------------------------------------------------------
    inh_mask = y_win <= (baseline_val - min_delta);

    cc = bwconncomp(inh_mask);
    valid_idx = [];

    for k = 1:cc.NumObjects
        idx = cc.PixelIdxList{k};
        if numel(idx) >= min_duration_bins
            valid_idx = [valid_idx; idx(:)];
        end
    end

    valid_idx = unique(valid_idx);
    response.n_valid_segments = sum(cellfun(@numel, cc.PixelIdxList) >= min_duration_bins);

    if isempty(valid_idx)
        return;
    end

    response.inh_present = true;

    % --------------------------------------------------------
    % Picco inibitorio = minimo all'interno dei segmenti validi
    % --------------------------------------------------------
    [peak_abs, loc] = min(y_win(valid_idx));
    true_idx = valid_idx(loc);

    response.peak_inh_abs = peak_abs;
    response.peak_inh_rel = peak_abs - baseline_val;   % valore negativo
    response.peak_inh_amp = baseline_val - peak_abs;   % valore positivo
    response.peak_inh_time = t_win(true_idx);

    % --------------------------------------------------------
    % Area thresholded: solo porzioni che superano il criterio
    % rispetto alla soglia baseline - min_delta
    % --------------------------------------------------------
    thr = baseline_val - min_delta;
    below_threshold = max(thr - y_win, 0);
    response.inh_area_thresholded = trapz(t_win, below_threshold);
end