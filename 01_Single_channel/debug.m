clearvars -except mean_baseline_common
close all
clc

% ============================================================
% FILE
% ============================================================
load('amplitude_component_analysis_global_baseline.mat', 'all_rows', ...
    'peak_window_mode', 'min_delta_from_baseline', 'min_duration_s')

filename = { ...
    '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat', ...
    '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat'};

cond_names = ["Gaze-on-center", "Gaze-on-target"];

% ============================================================
% PARAMETRI
% ============================================================
sets_plot = [1, 2, 3, 4, 5, 6];
n_sets = numel(sets_plot);

n_channels = 96;
bin_size = 0.02;
w = 15;
min_duration_bins = ceil(min_duration_s / bin_size);

% se true, mostra solo combinazioni con almeno una componente presente
show_only_present_components = false;

% ordine righe
all_rows = sortrows(all_rows, {'target_id','channel_global'});

% ============================================================
% CARICA CONDIZIONI
% ============================================================
cond_data = struct();

for d = 1:2
    S = load(filename{d});
    data = S.data;

    TS = data(1).Data(1).Interp(1).Task_states;
    state_names = string(TS(:,1));
    state_dur_s = cellfun(@(x) size(x,1) * bin_size, TS(:,2));
    state_onset_s = [0; cumsum(state_dur_s(1:end-1))];
    get_onset = @(name) state_onset_s(find(strcmpi(state_names, name), 1, 'first'));

    thisfile = string(filename{d});

    if contains(thisfile, "motor", "IgnoreCase", true)
        task_name = "Gaze-on-center";
        lines_t   = [get_onset("pres12"); get_onset("reach")];
        lines_lab = ["Target cue"; "Go cue"];

    elseif contains(thisfile, "controlled", "IgnoreCase", true)
        task_name = "Gaze-on-target";
        lines_t   = [get_onset("pres12"); get_onset("gaze"); get_onset("reach")];
        lines_lab = ["Target cue"; "Go cue - gaze"; "Go cue"];

    else
        error('File condizione non riconosciuto.');
    end

    cond_data(d).data = data;
    cond_data(d).task_name = task_name;
    cond_data(d).lines_t = lines_t;
    cond_data(d).lines_lab = lines_lab;
end

% ============================================================
% LOOP PRINCIPALE
% ============================================================
component_list = {'exc','inh'};
n_total_rows = height(all_rows);

for rr = 1:n_total_rows

    row = all_rows(rr,:);

    selected_target = row.target_id;
    selected_channel_global = row.channel_global;
    array = row.array_id;
    channel = row.channel_local;
    baseline_val = row.baseline;

    % --------------------------------------------------------
    % CALCOLA PSTH E PICCHI NELLE DUE CONDIZIONI
    % --------------------------------------------------------
    resp = repmat(empty_response_struct(), 1, 2);

    for d = 1:2
        resp(d) = compute_condition_response( ...
            cond_data(d).data, ...
            cond_data(d).lines_t, cond_data(d).lines_lab, ...
            array, channel, selected_target, ...
            sets_plot, n_sets, bin_size, w, ...
            peak_window_mode, baseline_val, ...
            min_delta_from_baseline, min_duration_bins);
    end

    % dati del target per scatter
    Tt = all_rows(all_rows.target_id == selected_target, :);

    for cc = 1:numel(component_list)

        scatter_type = component_list{cc};

        switch lower(scatter_type)
            case 'exc'
                Tplot = Tt(Tt.exc_any == true, :);
                x_all = Tplot.exc_center;
                y_all = Tplot.exc_target;

                x_sel = row.exc_center;
                y_sel = row.exc_target;

                x_present = row.exc_center_present;
                y_present = row.exc_target_present;

                scatter_title = sprintf('Target %d | Exc component', selected_target);
                xlab = 'Gaze-on-center exc';
                ylab = 'Gaze-on-target exc';

                if show_only_present_components && ~(x_present || y_present)
                    continue
                end

            case 'inh'
                Tplot = Tt(Tt.inh_any == true, :);
                x_all = Tplot.inh_center_amp;
                y_all = Tplot.inh_target_amp;

                x_sel = row.inh_center_amp;
                y_sel = row.inh_target_amp;

                x_present = row.inh_center_present;
                y_present = row.inh_target_present;

                scatter_title = sprintf('Target %d | Inh component', selected_target);
                xlab = 'Gaze-on-center inh amplitude';
                ylab = 'Gaze-on-target inh amplitude';

                if show_only_present_components && ~(x_present || y_present)
                    continue
                end

            otherwise
                error('scatter_type non riconosciuto.');
        end

        fprintf('\n========================================\n');
        fprintf('DEBUG COMPONENTI\n');
        fprintf('Progress: row %d/%d | component %s\n', rr, n_total_rows, upper(scatter_type));
        fprintf('Target: %d\n', selected_target);
        fprintf('Canale globale: %d\n', selected_channel_global);
        fprintf('Array: %d\n', array);
        fprintf('Canale locale: %d\n', channel);
        fprintf('Baseline globale: %.3f Hz\n', baseline_val);
        fprintf('Scatter point: (%.3f, %.3f)\n', x_sel, y_sel);
        fprintf('Presenza componente: center=%d | target=%d\n', x_present, y_present);
        fprintf('========================================\n');

        % ====================================================
        % FIGURA
        % ====================================================
        close all
        scr = get(groot, 'ScreenSize');
        figure('Color', 'w', 'Position', [60 60 min(1350,scr(3)-100) 850]);

        tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

        % ----------------------------------------------------
        % 1) Scatter con punto evidenziato
        % ----------------------------------------------------
        nexttile;
        hold on

        if ~isempty(Tplot)
            scatter(x_all, y_all, 40, 'filled', ...
                'MarkerFaceColor', [0.2 0.45 0.8], ...
                'MarkerFaceAlpha', 0.7);

            lims = [min([x_all; y_all]), max([x_all; y_all])];
            if isempty(lims) || any(isnan(lims))
                lims = [0 1];
            end
            if lims(1) == lims(2)
                lims = lims + [-1 1];
            end
            plot(lims, lims, '--k', 'LineWidth', 1.2);

            if ~isnan(x_sel) && ~isnan(y_sel)
                scatter(x_sel, y_sel, 140, 'o', ...
                    'MarkerFaceColor', [0.85 0.1 0.1], ...
                    'MarkerEdgeColor', 'k', ...
                    'LineWidth', 1.5);

                text(x_sel, y_sel, sprintf('  ch %d', selected_channel_global), ...
                    'FontWeight', 'bold', ...
                    'VerticalAlignment', 'bottom');
            end

            xlim(lims);
            ylim(lims);

            txt_sc = sprintf('center present = %d\ntarget present = %d', x_present, y_present);
            text(lims(1) + 0.05*range(lims), lims(2) - 0.08*range(lims), txt_sc, ...
                'FontSize', 10, ...
                'BackgroundColor', 'w', ...
                'EdgeColor', 'k');
        else
            text(0.5, 0.5, 'Nessun punto nel pool scatter per questo target', ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
        end

        xlabel(xlab);
        ylabel(ylab);
        title(scatter_title);
        grid on
        axis square
        hold off

        % ----------------------------------------------------
        % 2) PSTH Gaze-on-center
        % ----------------------------------------------------
        nexttile;
        plot_condition_panel(resp(1), cond_data(1).task_name, baseline_val, min_delta_from_baseline);

        % ----------------------------------------------------
        % 3) PSTH Gaze-on-target
        % ----------------------------------------------------
        nexttile;
        plot_condition_panel(resp(2), cond_data(2).task_name, baseline_val, min_delta_from_baseline);

        % ----------------------------------------------------
        % 4) Testo riassuntivo
        % ----------------------------------------------------
        nexttile;
        axis off

        txt = {
            sprintf('Target = %d', selected_target)
            sprintf('Canale globale = %d', selected_channel_global)
            sprintf('Array = %d | Canale locale = %d', array, channel)
            ' '
            sprintf('Baseline globale = %.3f Hz', baseline_val)
            ' '
            sprintf('Center response_type = %s', string(resp(1).response_type))
            sprintf('Center exc peak = %.3f Hz @ %.3f s', resp(1).peak_exc_rel, resp(1).peak_exc_time)
            sprintf('Center inh peak = %.3f Hz @ %.3f s', resp(1).peak_inh_rel, resp(1).peak_inh_time)
            sprintf('Center exc used in scatter = %.3f', row.exc_center)
            sprintf('Center inh used in scatter = %.3f', row.inh_center_amp)
            sprintf('Center exc present = %d | inh present = %d', row.exc_center_present, row.inh_center_present)
            ' '
            sprintf('Target response_type = %s', string(resp(2).response_type))
            sprintf('Target exc peak = %.3f Hz @ %.3f s', resp(2).peak_exc_rel, resp(2).peak_exc_time)
            sprintf('Target inh peak = %.3f Hz @ %.3f s', resp(2).peak_inh_rel, resp(2).peak_inh_time)
            sprintf('Target exc used in scatter = %.3f', row.exc_target)
            sprintf('Target inh used in scatter = %.3f', row.inh_target_amp)
            sprintf('Target exc present = %d | inh present = %d', row.exc_target_present, row.inh_target_present)
            ' '
            sprintf('Scatter type = %s', upper(scatter_type))
            sprintf('Selected point = (%.3f, %.3f)', x_sel, y_sel)
            };

        text(0.02, 0.98, txt, 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'FontSize', 10.5);

        sgtitle(sprintf('Debug componenti | row %d/%d | target %d | ch %d | %s', ...
            rr, n_total_rows, selected_target, selected_channel_global, upper(scatter_type)), ...
            'FontWeight', 'bold');

        % ----------------------------------------------------
        % CONSOLE
        % ----------------------------------------------------
        fprintf('\n--- GAZE-ON-CENTER ---\n');
        fprintf('response_type = %s\n', resp(1).response_type);
        fprintf('peak_exc_rel = %.3f | peak_exc_time = %.3f\n', resp(1).peak_exc_rel, resp(1).peak_exc_time);
        fprintf('peak_inh_rel = %.3f | peak_inh_time = %.3f\n', resp(1).peak_inh_rel, resp(1).peak_inh_time);
        fprintf('used exc = %.3f | used inh = %.3f\n', row.exc_center, row.inh_center_amp);

        fprintf('\n--- GAZE-ON-TARGET ---\n');
        fprintf('response_type = %s\n', resp(2).response_type);
        fprintf('peak_exc_rel = %.3f | peak_exc_time = %.3f\n', resp(2).peak_exc_rel, resp(2).peak_exc_time);
        fprintf('peak_inh_rel = %.3f | peak_inh_time = %.3f\n', resp(2).peak_inh_rel, resp(2).peak_inh_time);
        fprintf('used exc = %.3f | used inh = %.3f\n', row.exc_target, row.inh_target_amp);

        % ----------------------------------------------------
        % AVANZAMENTO MANUALE
        % ----------------------------------------------------
        str = input(sprintf(['Plot completato.\n' ...
                             'Premi INVIO per continuare, oppure digita q e premi INVIO per uscire: ']), 's');

        if strcmpi(str, 'q')
            disp('Interruzione richiesta dall''utente.');
            return;
        end
    end
end

disp('Fine debug di tutte le combinazioni.');

% ============================================================
% FUNZIONI LOCALI
% ============================================================
function response = compute_condition_response( ...
    data_cond, lines_t, lines_lab, ...
    array, channel, target_id, ...
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
            firing_rate = [firing_rate, M_spikes ./ bin_size]; %#ok<AGROW>
        end
    end

    if isempty(firing_rate)
        response = empty_response_struct();
        return;
    end

    mean_fr = mean(firing_rate, 2);
    psth = smoothdata(mean_fr, 'gaussian', w);
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

    response.t = t;
    response.psth = psth;
    response.mean_fr = mean_fr;
    response.go_time = go_time;
    response.win_start = win_start;
    response.win_end = win_end;
    response.lines_t = lines_t;
    response.lines_lab = lines_lab;
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
    out.t = [];
    out.psth = [];
    out.mean_fr = [];
    out.go_time = NaN;
    out.win_start = NaN;
    out.win_end = NaN;
    out.lines_t = [];
    out.lines_lab = strings(0,1);
end

function plot_condition_panel(resp, title_txt, baseline_val, thr)
    hold on

    if isempty(resp.t)
        text(0.5, 0.5, 'Nessun dato', ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center');
        title(char(title_txt));
        axis off
        return
    end

    plot(resp.t, resp.mean_fr, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
    plot(resp.t, resp.psth, 'k', 'LineWidth', 1.8);

    yline(baseline_val, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);
    yline(baseline_val + thr, ':', 'Color', [0.5 0.5 0.5]);
    yline(baseline_val - thr, ':', 'Color', [0.5 0.5 0.5]);

    yl = ylim;
    patch([resp.win_start resp.win_end resp.win_end resp.win_start], ...
          [yl(1) yl(1) yl(2) yl(2)], ...
          [0.85 0.9 1], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

    plot(resp.t, resp.mean_fr, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
    plot(resp.t, resp.psth, 'k', 'LineWidth', 1.8);

    for k = 1:numel(resp.lines_t)
        xline(resp.lines_t(k), '--k', 'HandleVisibility', 'off');
    end

    if ~isnan(resp.peak_exc_abs)
        scatter(resp.peak_exc_time, resp.peak_exc_abs, 80, 'r', 'filled');
    end
    if ~isnan(resp.peak_inh_abs)
        scatter(resp.peak_inh_time, resp.peak_inh_abs, 80, 'b', 'filled');
    end

    xlabel('Time (s)');
    ylabel('Firing rate (Hz)');
    title(char(title_txt));
    grid on
    hold off
end