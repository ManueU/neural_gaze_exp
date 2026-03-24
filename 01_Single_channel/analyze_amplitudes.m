clearvars -except mean_baseline_common
close all
clc

% ============================================================
% FILE
% ============================================================
load('responsive_channels_free_gaze.mat', 'responsive_channels', 'results')

filename = { ...
    '../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat', ...
    '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat', ...
    '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat'};

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

% 'planning'  = Go cue - 0.5 s --> Go cue
% 'reaching'  = Go cue         --> Go cue + 3.0 s
% 'full'      = Go cue - 0.5 s --> Go cue + 3.0 s
peak_window_mode = 'full';

min_delta_from_baseline = 12;   % Hz
min_duration_s = 0.06;          % s
min_duration_bins = ceil(min_duration_s / bin_size);

show_summary_figures = true;

if ~exist('mean_baseline_common','var') || isempty(mean_baseline_common)
    error('Serve mean_baseline_common nel workspace.');
end

% ============================================================
% CARICAMENTO CONDIZIONI
% ============================================================
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

        % ----------------------------------------------------
        % BASELINE GLOBALE
        % ----------------------------------------------------
        if ch_global > numel(mean_baseline_common)
            continue;
        end

        baseline_val = mean_baseline_common(ch_global);
        if isnan(baseline_val)
            continue;
        end

        % ----------------------------------------------------
        % FREE-GAZE: recupero risultati già salvati
        % ----------------------------------------------------
        idx_free = find(results.target_id == target_id & ...
                        results.channel_global == ch_global, 1, 'first');

        if isempty(idx_free)
            continue;
        end

        free_class_raw = string(results.response_type(idx_free));
        free_exc = results.peak_exc_rel(idx_free);
        free_inh = results.peak_inh_rel(idx_free);

        % ----------------------------------------------------
        % GAZE-ON-CENTER
        % ----------------------------------------------------
        resp_center = compute_condition_response( ...
            cond_data(2).data, ...
            cond_data(2).lines_t, cond_data(2).lines_lab, ...
            array, channel, ch_global, target_id, ...
            sets_plot, n_sets, bin_size, w, ...
            peak_window_mode, baseline_val, ...
            min_delta_from_baseline, min_duration_bins);

        % ----------------------------------------------------
        % GAZE-ON-TARGET
        % ----------------------------------------------------
        resp_target = compute_condition_response( ...
            cond_data(3).data, ...
            cond_data(3).lines_t, cond_data(3).lines_lab, ...
            array, channel, ch_global, target_id, ...
            sets_plot, n_sets, bin_size, w, ...
            peak_window_mode, baseline_val, ...
            min_delta_from_baseline, min_duration_bins);

        % ----------------------------------------------------
        % COMPONENTE ECCITATORIA
        % Se assente in una condizione -> 0
        % ----------------------------------------------------
        exc_center_present = ~isnan(resp_center.peak_exc_rel);
        exc_target_present = ~isnan(resp_target.peak_exc_rel);

        if exc_center_present
            exc_center = resp_center.peak_exc_rel;
        else
            exc_center = 0;
        end

        if exc_target_present
            exc_target = resp_target.peak_exc_rel;
        else
            exc_target = 0;
        end

        exc_any = exc_center_present || exc_target_present;
        exc_both_present = exc_center_present && exc_target_present;

        % ----------------------------------------------------
        % COMPONENTE INIBITORIA
        % Se assente in una condizione -> 0
        % Uso ampiezza positiva
        % ----------------------------------------------------
        inh_center_present = ~isnan(resp_center.peak_inh_rel);
        inh_target_present = ~isnan(resp_target.peak_inh_rel);

        if inh_center_present
            inh_center_amp = abs(resp_center.peak_inh_rel);
        else
            inh_center_amp = 0;
        end

        if inh_target_present
            inh_target_amp = abs(resp_target.peak_inh_rel);
        else
            inh_target_amp = 0;
        end

        inh_any = inh_center_present || inh_target_present;
        inh_both_present = inh_center_present && inh_target_present;

        % ----------------------------------------------------
        % SALVA
        % ----------------------------------------------------
        new_row = table( ...
            target_id, ch_global, array, channel, baseline_val, ...
            free_class_raw, free_exc, free_inh, ...
            string(resp_center.response_type), ...
            resp_center.peak_exc_rel, resp_center.peak_exc_time, ...
            resp_center.peak_inh_rel, resp_center.peak_inh_time, ...
            string(resp_target.response_type), ...
            resp_target.peak_exc_rel, resp_target.peak_exc_time, ...
            resp_target.peak_inh_rel, resp_target.peak_inh_time, ...
            exc_center_present, exc_target_present, exc_both_present, exc_any, ...
            exc_center, exc_target, ...
            inh_center_present, inh_target_present, inh_both_present, inh_any, ...
            inh_center_amp, inh_target_amp, ...
            'VariableNames', { ...
            'target_id', 'channel_global', 'array_id', 'channel_local', 'baseline', ...
            'free_class_raw', 'free_exc_rel', 'free_inh_rel', ...
            'center_class_raw', ...
            'center_exc_rel', 'center_exc_time', ...
            'center_inh_rel', 'center_inh_time', ...
            'target_class_raw', ...
            'target_exc_rel', 'target_exc_time', ...
            'target_inh_rel', 'target_inh_time', ...
            'exc_center_present', 'exc_target_present', 'exc_both_present', 'exc_any', ...
            'exc_center', 'exc_target', ...
            'inh_center_present', 'inh_target_present', 'inh_both_present', 'inh_any', ...
            'inh_center_amp', 'inh_target_amp'});

        all_rows = [all_rows; new_row]; %#ok<AGROW>
    end

    % ========================================================
    % SUMMARY PER TARGET
    % ========================================================
    Tt = all_rows(all_rows.target_id == target_id, :);

    if isempty(Tt)
        continue;
    end

    fprintf('\nTarget %d\n', target_id);
    fprintf('Canali con componente eccitatoria in almeno una condizione: %d\n', sum(Tt.exc_any));
    fprintf('Canali con componente eccitatoria in entrambe: %d\n', sum(Tt.exc_both_present));
    fprintf('Canali con componente inibitoria in almeno una condizione: %d\n', sum(Tt.inh_any));
    fprintf('Canali con componente inibitoria in entrambe: %d\n', sum(Tt.inh_both_present));

    % ========================================================
    % STATISTICA CORRETTA:
    % solo casi in cui la componente è presente in entrambe
    % ========================================================
    Te_stats = Tt(Tt.exc_center_present == 1 & Tt.exc_target_present == 1, :);
    Ti_stats = Tt(Tt.inh_center_present == 1 & Tt.inh_target_present == 1, :);

    if height(Te_stats) >= 3
        p_exc = signrank(Te_stats.exc_center, Te_stats.exc_target);
        mean_diff_exc = mean(Te_stats.exc_target - Te_stats.exc_center);
        fprintf('Target %d | signrank exc: p = %.4g | N = %d | mean diff = %.3f Hz\n', ...
            target_id, p_exc, height(Te_stats), mean_diff_exc);
    else
        p_exc = NaN;
        fprintf('Target %d | exc: N insufficiente\n', target_id);
    end

    if height(Ti_stats) >= 3
        p_inh = signrank(Ti_stats.inh_center_amp, Ti_stats.inh_target_amp);
        mean_diff_inh = mean(Ti_stats.inh_target_amp - Ti_stats.inh_center_amp);
        fprintf('Target %d | signrank inh: p = %.4g | N = %d | mean diff = %.3f Hz\n', ...
            target_id, p_inh, height(Ti_stats), mean_diff_inh);
    else
        p_inh = NaN;
        fprintf('Target %d | inh: N insufficiente\n', target_id);
    end

    % --------------------------------------------------------
    % FIGURA RIASSUNTIVA
    % Scatter: tutti i casi con componente presente in almeno una
    % Statistica mostrata: solo casi presenti in entrambe
    % --------------------------------------------------------
    if show_summary_figures
        scr = get(groot, 'ScreenSize');
        figure('Color', 'w', 'Position', [80 80 min(1200,scr(3)-100) 500]);
        tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');

        % =====================================================
        % ---- 1) scatter eccitazione + statistica
        % =====================================================
        nexttile; hold on
        Te_plot = Tt(Tt.exc_any == true, :);

        if ~isempty(Te_plot)
            scatter(Te_plot.exc_center, Te_plot.exc_target, 40, 'filled');

            lims = [min([Te_plot.exc_center; Te_plot.exc_target]), max([Te_plot.exc_center; Te_plot.exc_target])];
            if lims(1) == lims(2)
                lims = lims + [-1 1];
            end

            plot(lims, lims, '--k');
            xlim(lims); ylim(lims);

            if height(Te_stats) >= 3
                txt = sprintf('p = %.3g\nN = %d\n(both present)', p_exc, height(Te_stats));
            else
                txt = sprintf('N = %d (no test)\n(both present)', height(Te_stats));
            end

            text(lims(1) + 0.05*range(lims), ...
                 lims(2) - 0.1*range(lims), ...
                 txt, ...
                 'FontSize', 11, ...
                 'BackgroundColor', 'w', ...
                 'EdgeColor', 'k');
        else
            text(0.5, 0.5, 'Nessuna componente eccitatoria', ...
                'Units','normalized', 'HorizontalAlignment','center');
        end

        xlabel('Gaze-on-center exc');
        ylabel('Gaze-on-target exc');
        title(sprintf('Target %d | Exc component', target_id));
        grid on; axis square; hold off

        % =====================================================
        % ---- 2) scatter inibizione + statistica
        % =====================================================
        nexttile; hold on
        Ti_plot = Tt(Tt.inh_any == true, :);

        if ~isempty(Ti_plot)
            scatter(Ti_plot.inh_center_amp, Ti_plot.inh_target_amp, 40, 'filled');

            lims = [min([Ti_plot.inh_center_amp; Ti_plot.inh_target_amp]), max([Ti_plot.inh_center_amp; Ti_plot.inh_target_amp])];
            if lims(1) == lims(2)
                lims = lims + [-1 1];
            end

            plot(lims, lims, '--k');
            xlim(lims); ylim(lims);

            if height(Ti_stats) >= 3
                txt = sprintf('p = %.3g\nN = %d\n(both present)', p_inh, height(Ti_stats));
            else
                txt = sprintf('N = %d (no test)\n(both present)', height(Ti_stats));
            end

            text(lims(1) + 0.05*range(lims), ...
                 lims(2) - 0.1*range(lims), ...
                 txt, ...
                 'FontSize', 11, ...
                 'BackgroundColor', 'w', ...
                 'EdgeColor', 'k');
        else
            text(0.5, 0.5, 'Nessuna componente inibitoria', ...
                'Units','normalized', 'HorizontalAlignment','center');
        end

        xlabel('Gaze-on-center inh amplitude');
        ylabel('Gaze-on-target inh amplitude');
        title(sprintf('Target %d | Inh component', target_id));
        grid on; axis square; hold off
    end
end

% ============================================================
% SALVATAGGIO
% ============================================================
save('amplitude_component_analysis_global_baseline.mat', 'all_rows', ...
    'peak_window_mode', 'min_delta_from_baseline', 'min_duration_s');

disp('File salvato: amplitude_component_analysis_global_baseline.mat');

% ============================================================
% SUMMARY GLOBALE
% ============================================================
fprintf('\n========================================\n');
fprintf('SUMMARY GLOBALE COMPONENTI\n');
fprintf('========================================\n');
fprintf('Canali-target con componente eccitatoria in almeno una condizione: %d\n', sum(all_rows.exc_any));
fprintf('Canali-target con componente inibitoria in almeno una condizione: %d\n', sum(all_rows.inh_any));

% ============================================================
% FUNZIONI LOCALI
% ============================================================
function response = compute_condition_response( ...
    data_cond, lines_t, lines_lab, ...
    array, channel, ch_global, target_id, ...
    sets_plot, n_sets, bin_size, w, ...
    peak_window_mode, baseline_val, ...
    min_delta_from_baseline, min_duration_bins)

    firing_rate = [];

    for set_ = 1:n_sets
        set_id = sets_plot(set_);
        idx = find(([data_cond(set_id).Data(array).Interp.Target_ID] == target_id) & ...
                   ([data_cond(set_id).Data(array).Interp.Excluded] == 0));

        for j = 1:length(idx)
            M_spikes = data_cond(set_id).Data(array).Interp(idx(j)).Trial(:, channel);
            firing_rate = [firing_rate, M_spikes ./ bin_size]; 
        end
    end

    if isempty(firing_rate)
        response = empty_response_struct();
        return;
    end

    psth = smoothdata(mean(firing_rate, 2), 'gaussian', w);
    t = (0:numel(psth)-1) * bin_size;

    idx_go = find(strcmp(lines_lab, "Go cue"), 1);
    if isempty(idx_go)
        response = empty_response_struct();
        return;
    end

    go_time = lines_t(idx_go);

    switch lower(peak_window_mode)
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
            error('peak_window_mode non riconosciuto.');
    end

    win_start = max(win_start, t(1));
    win_end   = min(win_end,   t(end));

    idx_win = t >= win_start & t <= win_end;
    if ~any(idx_win)
        response = empty_response_struct();
        return;
    end

    t_win = t(idx_win);
    y_win = psth(idx_win);

    response = find_response_peaks( ...
        t_win, y_win, baseline_val, ...
        min_delta_from_baseline, min_duration_bins);
end

function out = empty_response_struct()
    out = struct();
    out.baseline = NaN;
    out.threshold = NaN;
    out.show_exc = false;
    out.show_inh = false;
    out.peak_exc_abs = NaN;
    out.peak_exc_rel = NaN;
    out.peak_exc_time = NaN;
    out.peak_inh_abs = NaN;
    out.peak_inh_rel = NaN;
    out.peak_inh_time = NaN;
    out.valid_exc = [];
    out.valid_inh = [];
    out.response_type = "none";
end