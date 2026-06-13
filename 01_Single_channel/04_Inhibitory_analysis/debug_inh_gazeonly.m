clearvars -except mean_baseline_common
close all
clc

% ============================================================
% INTERACTIVE DEBUG VIEWER
% gaze-only inhibitory peak detection
% Scorre tutti i channel-target selezionati in free-gaze
% ============================================================

% ============================================================
% FILE
% ============================================================
filename_gaze_only = '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat';
filename_free_gaze = 'responsive_channels_free_gaze.mat';

% ============================================================
% PARAMETRI
% ============================================================
sets_plot = [1, 2, 3, 4, 5, 6];
n_sets = numel(sets_plot);

n_arrays   = 2;
n_channels = 96;
bin_size   = 0.02;
w          = 15;

win_start_rel = -1.0;
win_end_rel   =  1.0;

min_delta_from_baseline = 12;   % Hz
min_duration_s = 0.06;          % s
min_duration_bins = ceil(min_duration_s / bin_size);

array_names = ["medial", "lateral"];

show_full_psth = true;

% ============================================================
% CHECK BASELINE
% ============================================================
if ~exist('mean_baseline_common','var') || isempty(mean_baseline_common)
    error('Serve mean_baseline_common nel workspace.');
end

% ============================================================
% LOAD FREE-GAZE RESPONSIVE CHANNELS
% ============================================================
load(filename_free_gaze, 'responsive_channels');

if ~exist('responsive_channels','var') || isempty(responsive_channels)
    error('responsive_channels non disponibile.');
end

% ============================================================
% COSTRUZIONE LISTA channel-target DA SCORRERE
% Una riga per ogni combinazione target-canale responsiva in free-gaze
% ============================================================
pairs = [];

for t_id = 1:numel(responsive_channels)
    target_id = responsive_channels(t_id).target_id;
    chs = responsive_channels(t_id).channels_global(:);

    if isempty(chs)
        continue;
    end

    pairs = [pairs; [repmat(target_id, numel(chs), 1), chs]]; %#ok<AGROW>
end

if isempty(pairs)
    error('Nessuna combinazione channel-target trovata.');
end

% ordina per target poi canale
pairs = sortrows(pairs, [1 2]);

fprintf('\n============================================\n');
fprintf('INTERACTIVE DEBUG VIEWER\n');
fprintf('============================================\n');
fprintf('Totale combinazioni channel-target: %d\n', size(pairs,1));
fprintf('Tasti:\n');
fprintf('  ->  / space / enter : prossimo\n');
fprintf('  <-  / b             : precedente\n');
fprintf('  q / ESC             : esci\n');

% ============================================================
% LOAD GAZE-ONLY DATA
% ============================================================
S = load(filename_gaze_only);
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

% ============================================================
% FIGURA
% ============================================================
hfig = figure('Color','w','Position',[60 60 1450 560], ...
              'Name','Gaze-only inhibitory debug viewer', ...
              'NumberTitle','off');

idx_pair = 1;

while ishandle(hfig)

    % --------------------------------------------------------
    % Estrai pair corrente
    % --------------------------------------------------------
    target_id = pairs(idx_pair, 1);
    ch_global = pairs(idx_pair, 2);

    array_id = ceil(ch_global / n_channels);
    channel_local = ch_global - (array_id - 1) * n_channels;

    if ch_global > numel(mean_baseline_common)
        baseline_val = NaN;
    else
        baseline_val = mean_baseline_common(ch_global);
    end

    % --------------------------------------------------------
    % Raccogli trial
    % --------------------------------------------------------
    firing_rate = [];
    trial_counter = 0;

    for set_ = 1:n_sets
        set_id = sets_plot(set_);
        idx_trials = find(([data(set_id).Data(array_id).Interp.Target_ID] == target_id) & ...
                          ([data(set_id).Data(array_id).Interp.Excluded] == 0));

        for j = 1:length(idx_trials)
            M_spikes = data(set_id).Data(array_id).Interp(idx_trials(j)).Trial(:, channel_local);
            firing_rate = [firing_rate, M_spikes ./ bin_size]; %#ok<AGROW>
            trial_counter = trial_counter + 1;
        end
    end

    clf(hfig);

    if isempty(firing_rate) || isnan(baseline_val)
        annotation('textbox', [0.2 0.4 0.6 0.2], ...
            'String', sprintf(['Pair %d / %d\nTarget %d | ch %d\n' ...
                               'Nessun trial o baseline non valida'], ...
                               idx_pair, size(pairs,1), target_id, ch_global), ...
            'FitBoxToText','on', 'HorizontalAlignment','center', ...
            'FontSize', 14, 'LineStyle','none');
        drawnow;
    else
        % ----------------------------------------------------
        % PSTH
        % ----------------------------------------------------
        psth = smoothdata(mean(firing_rate, 2), 'gaussian', w);
        t = (0:numel(psth)-1) * bin_size;

        this_win_start = gaze_time + win_start_rel;
        this_win_end   = gaze_time + win_end_rel;

        this_win_start = max(this_win_start, t(1));
        this_win_end   = min(this_win_end,   t(end));

        idx_win = t >= this_win_start & t <= this_win_end;

        if ~any(idx_win)
            annotation('textbox', [0.2 0.4 0.6 0.2], ...
                'String', sprintf(['Pair %d / %d\nTarget %d | ch %d\n' ...
                                   'Finestra di analisi vuota'], ...
                                   idx_pair, size(pairs,1), target_id, ch_global), ...
                'FitBoxToText','on', 'HorizontalAlignment','center', ...
                'FontSize', 14, 'LineStyle','none');
            drawnow;
        else
            t_win_abs = t(idx_win);
            y_win = psth(idx_win);
            t_win_rel = t_win_abs - gaze_time;

            dbg = debug_find_inhibitory_response( ...
                t_win_rel, y_win, baseline_val, ...
                min_delta_from_baseline, min_duration_bins);

            % ------------------------------------------------
            % LAYOUT
            % ------------------------------------------------
            if show_full_psth
                tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
            else
                tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
            end

            % --------------------------------------------
            % PANNELLO 1: PSTH intero
            % --------------------------------------------
            if show_full_psth
                nexttile; hold on

                yl_pad = 0.1 * range([psth(:); baseline_val]);
                if yl_pad == 0, yl_pad = 5; end

                yl0 = [min([psth(:); baseline_val]) - yl_pad, ...
                       max([psth(:); baseline_val]) + yl_pad];

                patch([this_win_start this_win_end this_win_end this_win_start], ...
                      [yl0(1) yl0(1) yl0(2) yl0(2)], ...
                      [0.92 0.92 0.92], 'EdgeColor','none', 'FaceAlpha',0.3);

                plot(t, psth, 'LineWidth', 1.6);
                xline(gaze_time, '--k', 'gaze', 'LabelVerticalAlignment','bottom');
                xline(this_win_start, ':k');
                xline(this_win_end, ':k');
                yline(baseline_val, '--', 'Baseline');
                yline(baseline_val - min_delta_from_baseline, '--', 'Threshold');

                ylim(yl0);
                xlabel('Time (s)');
                ylabel('Firing rate (Hz)');
                title(sprintf('Full PSTH | %s | ch %d (glob %d) | target %d', ...
                    array_names(array_id), channel_local, ch_global, target_id));
                grid on
                box off
                hold off
            end

            % --------------------------------------------
            % PANNELLO 2: finestra di detection
            % --------------------------------------------
            if show_full_psth
                nexttile;
            end
            hold on

            plot(t_win_rel, y_win, 'LineWidth', 1.8);
            yline(baseline_val, '--', 'Baseline', 'LineWidth', 1.1);
            yline(baseline_val - min_delta_from_baseline, '--', 'Threshold', 'LineWidth', 1.1);
            xline(0, '--k', 'gaze', 'LabelVerticalAlignment','bottom');

            % punti sotto threshold
            if any(dbg.inh_mask)
                scatter(t_win_rel(dbg.inh_mask), y_win(dbg.inh_mask), ...
                    22, 'filled', 'MarkerFaceAlpha', 0.35);
            end

            % segmenti validi
            for k = 1:dbg.n_valid_segments
                idxk = dbg.valid_segments{k};
                plot(t_win_rel(idxk), y_win(idxk), 'LineWidth', 3);
            end

            % picco selezionato
            if dbg.inh_present
                plot(dbg.peak_inh_time, dbg.peak_inh_abs, 'rv', ...
                    'MarkerFaceColor','r', 'MarkerSize',9);
            end

            xr = [win_start_rel, win_end_rel];
            xlim(xr);

            y_candidates = [y_win(:); baseline_val; baseline_val - min_delta_from_baseline];
            yr_pad = 0.1 * range(y_candidates);
            if yr_pad == 0, yr_pad = 5; end
            ylim([min(y_candidates)-yr_pad, max(y_candidates)+yr_pad]);

            if dbg.inh_present
                txt = sprintf([ ...
                    'Pair %d/%d\n' ...
                    'Trials = %d\n' ...
                    'Segments total = %d\n' ...
                    'Segments valid = %d\n' ...
                    'Peak abs = %.2f Hz\n' ...
                    'Peak amp = %.2f Hz\n' ...
                    'Peak time = %.2f s\n' ...
                    'Area total = %.2f Hz*s'], ...
                    idx_pair, size(pairs,1), trial_counter, ...
                    dbg.n_segments_total, dbg.n_valid_segments, ...
                    dbg.peak_inh_abs, dbg.peak_inh_amp, ...
                    dbg.peak_inh_time, dbg.inh_area_total);
            else
                txt = sprintf([ ...
                    'Pair %d/%d\n' ...
                    'Trials = %d\n' ...
                    'Segments total = %d\n' ...
                    'Segments valid = %d\n' ...
                    'NO VALID INHIBITION'], ...
                    idx_pair, size(pairs,1), trial_counter, ...
                    dbg.n_segments_total, dbg.n_valid_segments);
            end

            yr = ylim;
            text(xr(1)+0.03*range(xr), yr(2)-0.05*range(yr), txt, ...
                'FontSize',10, 'BackgroundColor','w', 'EdgeColor','k', ...
                'VerticalAlignment','top');

            xlabel('Time relative to gaze onset (s)');
            ylabel('Firing rate (Hz)');
            title('Detection window and selected inhibitory peak');
            grid on
            box off
            hold off

            drawnow;
        end
    end

    % --------------------------------------------------------
    % Aspetta input da tastiera
    % --------------------------------------------------------
    waitforbuttonpress;
    key = get(hfig, 'CurrentKey');

    switch lower(key)
        case {'rightarrow', 'space', 'return'}
            idx_pair = min(idx_pair + 1, size(pairs,1));

        case {'leftarrow', 'b'}
            idx_pair = max(idx_pair - 1, 1);

        case {'q', 'escape'}
            break

        otherwise
            % qualsiasi altro tasto: vai avanti
            idx_pair = min(idx_pair + 1, size(pairs,1));
    end
end

if ishandle(hfig)
    disp('Viewer chiuso.');
else
    disp('Figura chiusa manualmente.');
end

% ============================================================
% FUNZIONE LOCALE
% ============================================================
function out = debug_find_inhibitory_response(t_win, y_win, baseline_val, min_delta, min_duration_bins)

    out = struct();

    out.inh_present = false;
    out.peak_inh_abs = NaN;
    out.peak_inh_rel = NaN;
    out.peak_inh_amp = NaN;
    out.peak_inh_time = NaN;
    out.inh_area_total = NaN;
    out.inh_area_thresholded = NaN;

    out.inh_mask = [];
    out.valid_idx = [];
    out.valid_segments = {};
    out.n_segments_total = 0;
    out.n_valid_segments = 0;

    % area totale sotto baseline
    below_baseline = max(baseline_val - y_win, 0);
    out.inh_area_total = trapz(t_win, below_baseline);

    % threshold
    thr = baseline_val - min_delta;
    inh_mask = y_win <= thr;
    out.inh_mask = inh_mask;

    cc = bwconncomp(inh_mask);
    out.n_segments_total = cc.NumObjects;

    valid_idx = [];
    valid_segments = {};

    for k = 1:cc.NumObjects
        idx = cc.PixelIdxList{k};
        if numel(idx) >= min_duration_bins
            valid_idx = [valid_idx; idx(:)]; %#ok<AGROW>
            valid_segments{end+1} = idx(:); %#ok<AGROW>
        end
    end

    valid_idx = unique(valid_idx);

    out.valid_idx = valid_idx;
    out.valid_segments = valid_segments;
    out.n_valid_segments = numel(valid_segments);

    if isempty(valid_idx)
        out.inh_area_thresholded = trapz(t_win, max(thr - y_win, 0));
        return;
    end

    out.inh_present = true;

    [peak_abs, loc] = min(y_win(valid_idx));
    true_idx = valid_idx(loc);

    out.peak_inh_abs = peak_abs;
    out.peak_inh_rel = peak_abs - baseline_val;
    out.peak_inh_amp = baseline_val - peak_abs;
    out.peak_inh_time = t_win(true_idx);

    below_threshold = max(thr - y_win, 0);
    out.inh_area_thresholded = trapz(t_win, below_threshold);
end