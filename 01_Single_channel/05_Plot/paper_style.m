% =========================================================
% FIGURE SEPARATE STILE PAPER
%
% 1) una figura per condizione:
%    - motor
%    - controlled
%    - gaze-only
%
%    In ogni figura:
%    - PSTH/hist sopra
%    - raster sotto
%
% 2) una quarta figura con polar plot
%
% OPZIONI:
% - scelta array e canale
% - scelta target mostrato nei raster/PSTH
% - flag: 'hist' oppure 'psth'
% =========================================================

clearvars -except mean_baseline_common
close all
clc

%% =========================
% PARAMETRI UTENTE
% ==========================
sets_plot    = [1 2 3 4 5 6];
bin_size     = 0.02;      % s
n_channels   = 96;

% selezione canale
array_sel    = 2;         % 1=medial, 2=lateral
channel_sel  = 70;        % canale locale
ch_global    = (array_sel-1)*n_channels + channel_sel;

array_names = ["medial", "lateral"];

% target
target_ids        = 1:8;
target_show       = 1;    % target visualizzato in raster + psth/hist
target_angles_deg = [0 45 90 135 180 225 270 315];

% modalità plot: 'hist' oppure 'psth'
plot_mode = 'psth';
smooth_w  = 15;            % smoothing per PSTH

% finestra relativa per polar plot (in secondi)
win_rel = [-0.5 0.5];

% file condizioni
cond(1).name = 'Gaze-on-center';
cond(1).file = '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat';

cond(2).name = 'Gaze-on-target';
cond(2).file = '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat';

cond(3).name = 'Gaze-only';
cond(3).file = '../00_Data_extraction/BCI02_Session_0924/gaze_BCI02_exclUpdated.mat';

% stili polar plot
cond_colors.motor      = [0 0 0];
cond_colors.controlled = [0 0 0];
cond_colors.gazeonly   = [0 0 0];

cond_lines.motor       = '-';
cond_lines.controlled  = '--';
cond_lines.gazeonly    = ':';

%% =========================
% POLAR DATA
% ==========================
polar_activity = nan(numel(cond), numel(target_ids));

%% =========================
% FIGURE DELLE 3 CONDIZIONI
% ==========================
for c = 1:numel(cond)

    S = load(cond(c).file);
    data = S.data;

    if cond(c).name == "Gaze-only"
        t_ref = get_onset(data, bin_size, 'reach');   
    else 
        t_ref = get_onset(data, bin_size, 'gaze');  
    end
    analysis_win_s = t_ref + win_rel;

    [event_times, event_labels] = get_condition_events(data, bin_size, cond(c).name);

    [trial_mat, raster_x, raster_y] = collect_trials_for_target( ...
        data, sets_plot, array_sel, channel_sel, target_show, bin_size);

    if isempty(trial_mat)
        warning('Nessun trial trovato per condizione %s, target %d.', cond(c).name, target_show);
        continue
    end

    % Figura
    fig = figure('Color','w','Position',[100 100 850 700]);

    tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    title_str = sprintf('%s | Channel %d | Target %d', cond(c).name, ch_global, target_show);
    sgtitle(title_str, 'FontWeight','bold', 'FontSize',16);

    % =========================
    % ASSE SOPRA: HIST / PSTH
    % =========================
    ax1 = nexttile(1);
    hold(ax1,'on')

    t = (0:size(trial_mat,1)-1) * bin_size;
    rate = mean(trial_mat, 2) ./ bin_size;
    centered_rate = rate - mean_baseline_common(ch_global);

    switch lower(plot_mode)
        case 'hist'
            bar(ax1, t, centered_rate, 1.0, 'FaceColor', 'none', 'EdgeColor', 'k');

        case 'psth'
            centered_rate_s = smoothdata(centered_rate, 'gaussian', smooth_w);
            plot(ax1, t, centered_rate_s, 'k', 'LineWidth', 1.8);

        otherwise
            error('plot_mode deve essere ''hist'' oppure ''psth''.')
    end

    % linee eventi
    for k = 1:numel(event_times)
        xline(ax1, event_times(k), '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.0, 'HandleVisibility','off');
    end

    ylabel(ax1, 'Firing rate (Hz)')
    xlim(ax1, [0 t(end)])
    box(ax1,'off')
    ax1.XTickLabel = [];   % niente etichette x sopra

    % =========================
    % ASSE SOTTO: RASTER
    % =========================
    ax2 = nexttile(2);
    hold(ax2,'on')

    % raster come lineette verticali
    for i = 1:numel(raster_x)
        line(ax2, [raster_x(i) raster_x(i)], [raster_y(i)-0.2, raster_y(i)+0.2], 'Color', 'k', 'LineWidth', 0.7);
    end

    % linee eventi + etichette
    for k = 1:numel(event_times)
        xline(ax2, event_times(k), '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.0, 'HandleVisibility','off');
    end

    xlabel(ax2, 'Time (s)')
    ylabel(ax2, 'Trials')
    xlim(ax2, [0 t(end)])
    ylim(ax2, [0.5 max(1, size(trial_mat,2)+0.5)])
    set(ax2, 'YDir', 'reverse')   % opzionale: primo trial in alto come nei raster classici
    box(ax2,'off')

    % etichette eventi sotto il raster
    yl2 = ylim(ax2);
    y_text = yl2(2) + 0.08 * range(yl2);
    for k = 1:numel(event_times)
        text(ax2, event_times(k), y_text, event_labels(k), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize', 10, ...
            'Clipping','off');
    end

    % =========================
    % DATI POLAR
    % =========================
    for tt = 1:numel(target_ids)
        trial_mat_t = collect_trials_matrix_only( ...
            data, sets_plot, array_sel, channel_sel, target_ids(tt));

        if isempty(trial_mat_t)
            polar_activity(c,tt) = nan;
            continue
        end

        tvec = (0:size(trial_mat_t,1)-1) * bin_size;
        idx_win = tvec >= analysis_win_s(1) & tvec < analysis_win_s(2);

        if any(idx_win)
            polar_activity(c,tt) = (mean(trial_mat_t(idx_win,:), 'all') / bin_size) - mean_baseline_common(ch_global);
        else
            polar_activity(c,tt) = nan;
        end
    end
end

%% =========================
% FIGURA 4: POLAR PLOT
% ==========================
fig4 = figure('Color','w','Position',[200 120 700 700]);
pax = polaraxes(fig4);
hold(pax,'on')

theta = deg2rad(target_angles_deg);
theta_closed = [theta theta(1)];

for c = 1:numel(cond)
    r = polar_activity(c,:);
    if all(isnan(r))
        continue
    end
    r_closed = [r r(1)];

    switch lower(cond(c).name)
        case 'gaze-on-center'
            ls = cond_lines.motor;
        case 'gaze-on-target'
            ls = cond_lines.controlled;
        case 'gaze-only'
            ls = cond_lines.gazeonly;
        otherwise
            ls = '-';
    end

    polarplot(pax, theta_closed, r_closed, ...
        'k', 'LineWidth', 2.0, 'LineStyle', ls, ...
        'DisplayName', cond(c).name);
end

pax.ThetaZeroLocation = 'right';
pax.ThetaDir = 'counterclockwise';
title(pax, sprintf('Polar tuning | Array %s | Ch %d (global %d)', ...
    array_names(array_sel), channel_sel, ch_global), ...
    'FontWeight','bold', 'FontSize',16)
legend(pax, 'Location', 'southoutside')

%% =========================================================
% FUNZIONI LOCALI
% =========================================================

function [event_times, event_labels] = get_condition_events(data, bin_size, cond_name)

    event_times = [];
    event_labels = strings(0,1);

    found = false;
    for s = 1:numel(data)
        for a = 1:numel(data(s).Data)
            if isfield(data(s).Data(a),'Interp') && ~isempty(data(s).Data(a).Interp)
                TS = data(s).Data(a).Interp(1).Task_states;
                found = true;
                break
            end
        end
        if found, break; end
    end

    if ~found
        warning('Task_states non trovati.');
        return
    end

    state_names = string(TS(:,1));
    state_dur_s = cellfun(@(x) size(x,1)*bin_size, TS(:,2));
    state_onset_s = [0; cumsum(state_dur_s(1:end-1))];

    get_onset = @(name) get_state_onset(state_names, state_onset_s, name);

    switch lower(cond_name)
        case 'motor'
            event_times  = [get_onset('pres12'), get_onset('reach')];
            event_labels = ["Target", "Go"];

        case 'controlled'
            event_times  = [get_onset('pres12'), get_onset('gaze'), get_onset('reach')];
            event_labels = ["Target", "Saccade", "Go"];

        case 'gaze-only'
            event_times  = [get_onset('pres12'), get_onset('gaze')];
            event_labels = ["Target", "Saccade"];

        otherwise
            event_times  = [get_onset('pres12')];
            event_labels = ["Target"];
    end

    valid = ~isnan(event_times);
    event_times = event_times(valid);
    event_labels = event_labels(valid);
end

function t = get_state_onset(state_names, state_onset_s, name)
    idx = find(strcmpi(state_names, name), 1, 'first');
    if isempty(idx)
        t = nan;
    else
        t = state_onset_s(idx);
    end
end

function [trial_mat, raster_x, raster_y] = collect_trials_for_target( ...
    data, sets_plot, array_sel, channel_sel, target_id, bin_size)

    trial_mat = [];
    raster_x = [];
    raster_y = [];

    row_id = 0;

    for iset = 1:numel(sets_plot)
        s = sets_plot(iset);
        if s > numel(data)
            continue
        end

        trials = data(s).Data(array_sel).Interp;
        idx = find([trials.Target_ID] == target_id & [trials.Excluded] == 0);

        for j = 1:numel(idx)
            spk = trials(idx(j)).Trial(:,channel_sel);

            if isempty(spk)
                continue
            end

            if isempty(trial_mat)
                trial_mat = spk;
            else
                L = min(size(trial_mat,1), numel(spk));
                trial_mat = trial_mat(1:L,:);
                spk = spk(1:L);
                trial_mat(:,end+1) = spk;
            end

            row_id = row_id + 1;

            st = find(spk > 0);
            if ~isempty(st)
                raster_x = [raster_x; (st-1)*bin_size];
                raster_y = [raster_y; row_id*ones(numel(st),1)]; 
            end
        end
    end
end

function trial_mat = collect_trials_matrix_only(data, sets_plot, array_sel, channel_sel, target_id)

    trial_mat = [];

    for iset = 1:numel(sets_plot)
        s = sets_plot(iset);
        if s > numel(data)
            continue
        end

        trials = data(s).Data(array_sel).Interp;
        idx = find([trials.Target_ID] == target_id & [trials.Excluded] == 0);

        for j = 1:numel(idx)
            spk = trials(idx(j)).Trial(:,channel_sel);

            if isempty(spk)
                continue
            end

            if isempty(trial_mat)
                trial_mat = spk;
            else
                L = min(size(trial_mat,1), numel(spk));
                trial_mat = trial_mat(1:L,:);
                spk = spk(1:L);
                trial_mat(:,end+1) = spk;
            end
        end
    end
end