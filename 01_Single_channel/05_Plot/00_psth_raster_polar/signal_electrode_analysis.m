%% Single-electrode directional tuning analysis
% Extracts target-specific trials, plots PSTH/raster and computes polar
% tuning for reach- and gaze-aligned epochs.

clearvars -except mean_baseline_common
close all
clc

%% Parameters
sets_plot = 1:6;
bin_size = 0.02;
n_channels = 96;

array_sel = 2;
channel_list = 70;
array_names = ["medial", "lateral"];

target_ids = 1:8;
target_angles_deg = [0 45 90 135 180 225 270 315];
target_colors = [
    0.80 0.20 0.20   % 0°   rosso
    0.90 0.55 0.15   % 45°  arancione
    0.85 0.80 0.20   % 90°  giallo
    0.35 0.70 0.25   % 135° verde
    0.10 0.65 0.55   % 180° turchese
    0.20 0.55 0.85   % 225° azzurro
    0.35 0.35 0.80   % 270° blu
    0.65 0.30 0.75   % 315° viola
];

smooth_w = 15;
win_rel_move = [-0.5 1.0];
win_rel_gaze = [-0.5 0.5];
display_mode = "byTarget";   % "allTargets" or "byTarget"

%% Conditions
cond = struct( ...
    'name', {'Free-gaze', 'Gaze-on-center', 'Gaze-on-target', 'Gaze-only'}, ...
    'code', {'free-gaze', 'motor', 'controlled', 'gaze-only'}, ...
    'file', { ...
        '../../../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat', ...
        '../../../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat', ...
        '../../../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat', ...
        '../../../00_Data_extraction/BCI02_Session_0924/gaze_BCI02_exclUpdated.mat'}, ...
    'line_style', {'-.', '-', '--', ':'});

reach_cond_idx = [1 2 3];
gaze_cond_idx  = [3 4];

%% Load data and events
fprintf('Loading data...\n');
t_load = tic;

for c = 1:numel(cond)
    S = load(cond(c).file);
    cond(c).data = S.data;
    [cond(c).event_times, cond(c).event_labels] = get_condition_events(cond(c).data, bin_size, cond(c).code);
end

fprintf('Loading completed in %.2f s\n\n', toc(t_load));

%% Main loop
for channel_sel = channel_list
    t_channel = tic;
    ch_global = (array_sel - 1) * n_channels + channel_sel;
    baseline_hz = mean_baseline_common(ch_global);

    fprintf('\nProcessing %s, channel %d, global %d...\n', array_names(array_sel), channel_sel, ch_global);

    trial_cache = build_trial_cache(cond, target_ids, sets_plot, array_sel, channel_sel, bin_size);

    polar_reach = compute_polar_activity(trial_cache, cond, reach_cond_idx, target_ids, ...
        bin_size, win_rel_move, baseline_hz, 'reach');
    polar_gaze = compute_polar_activity(trial_cache, cond, gaze_cond_idx, target_ids, ...
        bin_size, win_rel_gaze, baseline_hz, 'gaze');

    plot_psth_rasters(display_mode, trial_cache, cond, target_ids, target_colors, ...
        channel_sel, bin_size, smooth_w, baseline_hz);

    plot_polar_summary(polar_reach, polar_gaze, cond, reach_cond_idx, gaze_cond_idx, ...
        target_angles_deg, channel_sel);

    fprintf('Channel %d completed in %.2f s\n', channel_sel, toc(t_channel));
    fprintf('Press ENTER for next channel.\n\n');

    pause
    close all
end

%% Local functions
function plot_psth_rasters(display_mode, trial_cache, cond, target_ids, target_colors, channel_sel, bin_size, smooth_w, baseline_hz)

    t_plot = tic;
    psth_ylim = compute_psth_ylim(trial_cache, bin_size, smooth_w, baseline_hz);

    switch display_mode
        case "allTargets"
            fig = figure('Color','w');
            tl = tiledlayout(fig, 2, numel(cond), 'TileSpacing','compact', 'Padding','compact');

            for c = 1:numel(cond)
                plot_condition_all_targets(nexttile(tl, c), nexttile(tl, numel(cond) + c), ...
                    trial_cache(c,:), target_ids, target_colors, bin_size, smooth_w, ...
                    baseline_hz, cond(c), psth_ylim, c);
            end

        case "byTarget"
            for tt = 1:numel(target_ids)
                fig = figure('Color','w');
                tl = tiledlayout(fig, 2, numel(cond), 'TileSpacing','compact', 'Padding','compact');
                color = target_colors(tt, :); 
                for c = 1:numel(cond)
                    cache = trial_cache{c,tt};
                    plot_condition_single_target(nexttile(tl, c), nexttile(tl, numel(cond) + c), ...
                        cache, bin_size, smooth_w, baseline_hz, cond(c), color, psth_ylim, c);
                end
            end

        otherwise
            error('display_mode must be "allTargets" or "byTarget".');
    end

    fprintf('PSTH/raster channel %d: %.2f s\n', channel_sel, toc(t_plot));
end

function plot_condition_all_targets(ax_psth, ax_raster, trial_cache_row, target_ids, ...
        target_colors, bin_size, smooth_w, baseline_hz, cond_info, psth_ylim, cond_idx)

    prepare_axes(ax_psth, ax_raster);
    trial_offset = 0;
    t_end = [];

    for tt = 1:numel(target_ids)
        cache = trial_cache_row{tt};
        if isempty(cache.trial_mat)
            continue
        end

        color = target_colors(tt,:);
        t = plot_psth(ax_psth, cache.trial_mat, bin_size, smooth_w, baseline_hz, color, ...
            sprintf('T%d', target_ids(tt)));
        plot_raster(ax_raster, cache.raster_x, cache.raster_y + trial_offset, color);

        trial_offset = trial_offset + size(cache.trial_mat, 2);
        t_end = t(end);
    end

    finish_condition_axes(ax_psth, ax_raster, cond_info, t_end, trial_offset, psth_ylim, cond_idx);
    % legend(ax_psth, 'Location','northwest', 'Box','off', 'FontSize',8);
end

function plot_condition_single_target(ax_psth, ax_raster, cache, bin_size, smooth_w, baseline_hz, cond_info, color, psth_ylim, cond_idx)
    prepare_axes(ax_psth, ax_raster);

    if isempty(cache.trial_mat)
        title(ax_psth, sprintf('%s | no trials', cond_info.name), 'Interpreter','none');
        axis(ax_psth,'off');
        axis(ax_raster,'off');
        return
    end

    t = plot_psth(ax_psth, cache.trial_mat, bin_size, smooth_w, baseline_hz, color, '');
    plot_raster(ax_raster, cache.raster_x, cache.raster_y, color);
    finish_condition_axes(ax_psth, ax_raster, cond_info, t(end), size(cache.trial_mat, 2), psth_ylim, cond_idx);
end

function prepare_axes(ax_psth, ax_raster)
    cla(ax_psth); cla(ax_raster);
    hold(ax_psth,'on'); hold(ax_raster,'on');
end

function t = plot_psth(ax, trial_mat, bin_size, smooth_w, baseline_hz, color, label)
    t = (0:size(trial_mat,1)-1) * bin_size;
    rate_hz = mean(trial_mat, 2) ./ bin_size - baseline_hz;
    rate_hz = smoothdata(rate_hz, 'gaussian', smooth_w);

    args = {'Color', color, 'LineWidth', 1.5};
    if strlength(string(label)) > 0
        args = [args, {'DisplayName', label}]; 
    end
    plot(ax, t, rate_hz, args{:});    
end

function plot_raster(ax, raster_x, raster_y, color)
    if isempty(raster_x)
        return
    end

    n_spikes = numel(raster_x);
    x = nan(3*n_spikes, 1);
    y = nan(3*n_spikes, 1);

    x(1:3:end) = raster_x;
    x(2:3:end) = raster_x;
    y(1:3:end) = raster_y - 0.25;
    y(2:3:end) = raster_y + 0.25;

    plot(ax, x, y, 'Color', color, 'LineWidth', 0.6, 'HandleVisibility','off');
end

function finish_condition_axes(ax_psth, ax_raster, cond_info, t_end, n_trials, psth_ylim, cond_idx)

    add_event_lines(ax_psth, cond_info.event_times);
    add_event_lines(ax_raster, cond_info.event_times);

    title(ax_psth, cond_info.name, 'Interpreter','none');
    if cond_idx == 1
        ylabel(ax_psth, 'FR - baseline (Hz)');
    else
        ylabel(ax_psth, '');
    end
    box(ax_psth,'off');
    ax_psth.XTickLabel = [];

    if ~isempty(psth_ylim)
        ylim(ax_psth, psth_ylim);
    end
    yline(ax_psth, 0, '-k');

    add_event_labels_between_axes(ax_psth, cond_info.event_times, cond_info.event_labels);

    xlabel(ax_raster, 'Time (s)');
    if cond_idx == 1
        ylabel(ax_raster, 'Trials');
    else
        ylabel(ax_raster, '');
    end
    set(ax_raster, 'YDir','reverse');
    box(ax_raster,'off');

    if ~isempty(t_end)
        xlim(ax_psth, [0 t_end]);
        xlim(ax_raster, [0 t_end]);
    end

    ylim(ax_raster, [0.5 max(1, n_trials + 0.5)]);
    yticklabels(ax_raster, []);
end

function psth_ylim = compute_psth_ylim(trial_cache, bin_size, smooth_w, baseline_hz)

    global_min = inf;
    global_max = -inf;

    for i = 1:numel(trial_cache)
        cache = trial_cache{i};

        if isempty(cache.trial_mat)
            continue
        end

        rate_hz = mean(cache.trial_mat, 2) ./ bin_size - baseline_hz;
        rate_hz = smoothdata(rate_hz, 'gaussian', smooth_w);

        global_min = min(global_min, min(rate_hz));
        global_max = max(global_max, max(rate_hz));
    end

    if isinf(global_min) || isinf(global_max)
        psth_ylim = [];
        return
    end

    margin = 0.08 * (global_max - global_min);

    if margin == 0
        margin = 1;
    end

    psth_ylim = [global_min - margin, global_max + margin];
end

function add_event_lines(ax, event_times)
    for k = 1:numel(event_times)
        xline(ax, event_times(k), '--', 'Color',[0.2 0.2 0.2], ...
            'LineWidth',1.0, 'HandleVisibility','off');
    end
end

function add_event_labels_between_axes(ax_psth, event_times, event_labels)

    yl = ylim(ax_psth);
    y_text = yl(1) - 0.08 * range(yl);

    for k = 1:numel(event_times)
        text(ax_psth, event_times(k), y_text, event_labels(k), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','top', ...
            'FontSize',9, ...
            'Clipping','off');
    end
end

function plot_polar_summary(polar_reach, polar_gaze, cond, reach_cond_idx, gaze_cond_idx, ...
        target_angles_deg, channel_sel)

    t_polar = tic;
    fig = figure('Color','w');

    pax1 = polaraxes('Parent', fig, 'Position',[0.07 0.18 0.38 0.70]);
    plot_polar_tuning(pax1, polar_reach, cond(reach_cond_idx), target_angles_deg);
    title(pax1, 'Reach-aligned', 'FontWeight','bold');

    pax2 = polaraxes('Parent', fig, 'Position',[0.57 0.18 0.38 0.70]);
    plot_polar_tuning(pax2, polar_gaze, cond(gaze_cond_idx), target_angles_deg);
    title(pax2, 'Gaze-aligned', 'FontWeight','bold');

    fprintf('Polar channel %d: %.2f s\n', channel_sel, toc(t_polar));
end

function plot_polar_tuning(pax, polar_activity, cond_subset, target_angles_deg)
    hold(pax,'on');

    theta = deg2rad(target_angles_deg);
    theta_closed = [theta theta(1)];

    for c = 1:numel(cond_subset)
        r = polar_activity(c,:);
        if all(isnan(r))
            continue
        end
        r = abs(r);
        polarplot(pax, theta_closed, [r r(1)], 'k', 'LineWidth',1.8, ...
            'LineStyle',cond_subset(c).line_style, 'DisplayName',cond_subset(c).name);
    end

    vals = abs(polar_activity(:));
    vals = vals(~isnan(vals));
    if ~isempty(vals)
        pax.RLim = [0 max(vals) + eps];
    end

    pax.ThetaZeroLocation = 'right';
    pax.ThetaDir = 'counterclockwise';
    legend(pax, 'Location','southoutside');
end
