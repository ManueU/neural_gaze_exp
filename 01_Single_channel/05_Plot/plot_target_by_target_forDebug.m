% clearvars
close all
clc

% ============================================================
% FILES
% ============================================================
data_file = '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat';

% ============================================================
% PARAMETRI
% ============================================================
sets_plot = [1, 2, 3, 4, 5, 6];

n_channels = 96;
bin_size = 0.02;
w = 15;

array_names = ["medial", "lateral"];
target_des = 1:8;

% finestra solo da visualizzare come riferimento
win_start_rel = -1.0;
win_end_rel   =  1.0;
show_window_lines = true;

% debug
debug_pause = true;              % premi un tasto tra un canale e il successivo
close_after_pause = true;        % chiude la figura dopo la pausa

% filtri canali
% channels_to_plot = [35
%     37
%     43
%     47
%     56
%     61
%     98
%    100
%    107
%    115
%    120
%    121
%    123
%    131
%    133
%    135
%    158
%    163
%    169
%    170
%    172
%    180
%    181
%    182
%    185
%    187];  
channels_to_plot = [182];
% esempio:
% channels_to_plot = [38 61 100 115 133 134 163];
% lascia [] per plottare tutti i canali

% se true usa full screen
maximize_figure = true;

% ============================================================
% LOAD DATA
% ============================================================
S = load(data_file);
if ~isfield(S, 'data')
    error('Nel file dati non trovo la variabile data.');
end
data = S.data;


% ============================================================
% LOAD PICCHI TROVATI IN PRECEDENZA
% ============================================================
peak_file = 'inhibitory_analysis_got.mat';

if exist(peak_file, 'file')
    P = load(peak_file, 'all_rows');
    has_peak_table = isfield(P, 'all_rows');
    if has_peak_table
        peak_table = P.all_rows;
    else
        peak_table = table();
    end
else
    warning('File %s non trovato: i picchi non verranno mostrati.', peak_file);
    has_peak_table = false;
    peak_table = table();
end


% ============================================================
% EVENT ONSET
% ============================================================
[~, baseName, ext] = fileparts(data_file);
ds_name = [baseName ext];

TS = data(1).Data(1).Interp(1).Task_states;
state_names = string(TS(:,1));
state_dur_s = cellfun(@(x) size(x,1) * bin_size, TS(:,2));
state_onset_s = [0; cumsum(state_dur_s(1:end-1))];
get_onset = @(name) state_onset_s(find(strcmpi(state_names, name), 1, 'first'));


if strcmp(ds_name, 'free-gaze_BCI02_exclUpdated.mat') || strcmp(ds_name, 'motor_BCI02_exclUpdated.mat') 
    time_event = get_onset("reach");
else 
    time_event = get_onset("gaze");
end 

% ============================================================
% CANALI DA PLOTTARE
% ============================================================
n_total_channels = numel(mean_baseline_common);

if isempty(channels_to_plot)
    channels_unique = 1:n_total_channels;
else
    channels_unique = channels_to_plot(:)';
    channels_unique = channels_unique(channels_unique >= 1 & channels_unique <= n_total_channels);
end

if isempty(channels_unique)
    error('Nessun canale valido da plottare.');
end

fprintf('\n============================================\n');
fprintf('VISUALIZZAZIONE SOLO SEGNALE + BASELINE\n');
fprintf('============================================\n');
fprintf('Canali da plottare: %d\n', numel(channels_unique));

% ============================================================
% LOOP PER CANALE
% ============================================================
for cc = 1:numel(channels_unique)

    ch_global = channels_unique(cc);
    array = ceil(ch_global / n_channels);
    channel = ch_global - (array - 1) * n_channels;

    if array < 1 || array > numel(array_names)
        fprintf('\n[SKIP] Canale globale %d non valido.\n', ch_global);
        continue;
    end

    baseline_val = mean_baseline_common(ch_global);
    if isnan(baseline_val)
        fprintf('\n[SKIP] Canale globale %d con baseline NaN.\n', ch_global);
        continue;
    end

    fprintf('\n============================================\n');
    fprintf('Canale %d / %d | %s | ch %d (glob %d)\n', ...
        cc, numel(channels_unique), array_names(array), channel, ch_global);
    fprintf('Baseline = %.2f Hz\n', baseline_val);
    fprintf('============================================\n');

    % --------------------------------------------------------
    % setup figura
    % --------------------------------------------------------
    if maximize_figure
        f = figure('Color','w','Units','normalized','OuterPosition',[0.02 0.05 0.96 0.9]);
    else
        f = figure('Color','w','Position',[100 100 1600 850]);
    end

    nT = numel(target_des);
    ncols = 4;
    nrows = ceil(nT / ncols);

    tl = tiledlayout(nrows, ncols, 'TileSpacing','compact', 'Padding','compact');
    title(tl, sprintf('%s | ch %d (glob %d) | baseline = %.2f Hz', ...
        array_names(array), channel, ch_global, baseline_val), ...
        'Interpreter', 'none');

    % --------------------------------------------------------
    % prima passata: ricostruisci tutti i target per asse y comune
    % --------------------------------------------------------
    cache = struct();
    YMIN_ALL = inf;
    YMAX_ALL = -inf;
    TMIN_ALL = inf;
    TMAX_ALL = -inf;

    for tt = 1:nT

        target_id = target_des(tt);

        cache(tt).target_id = target_id;
        cache(tt).has_data = false;
        cache(tt).t_rel = [];
        cache(tt).psth = [];

        fr_cells = {};
        idx_cell = 1;

        for set_id = sets_plot
            idx_trials = find(([data(set_id).Data(array).Interp.Target_ID] == target_id) & ...
                              ([data(set_id).Data(array).Interp.Excluded] == 0));

            for j = 1:numel(idx_trials)
                M_spikes = data(set_id).Data(array).Interp(idx_trials(j)).Trial(:, channel);
                fr_cells{idx_cell} = M_spikes ./ bin_size; 
                idx_cell = idx_cell + 1;
            end
        end

        if isempty(fr_cells)
            continue;
        end

        firing_rate = [fr_cells{:}];
        psth = smoothdata(mean(firing_rate, 2), 'gaussian', w);

        t = (0:numel(psth)-1) * bin_size;
        t_rel = t - time_event;

        cache(tt).has_data = true;
        cache(tt).t_rel = t_rel;
        cache(tt).psth = psth;

        YMIN_ALL = min(YMIN_ALL, min([psth(:); baseline_val]));
        YMAX_ALL = max(YMAX_ALL, max([psth(:); baseline_val]));
        TMIN_ALL = min(TMIN_ALL, t_rel(1));
        TMAX_ALL = max(TMAX_ALL, t_rel(end));
    end

    if isinf(YMIN_ALL) || isinf(YMAX_ALL)
        YMIN_ALL = baseline_val - 5;
        YMAX_ALL = baseline_val + 5;
    end

    if isinf(TMIN_ALL) || isinf(TMAX_ALL)
        TMIN_ALL = -1;
        TMAX_ALL = 1;
    end

    YMIN_ALL = YMIN_ALL - 2;
    YMAX_ALL = YMAX_ALL + 2;

    % --------------------------------------------------------
    % seconda passata: plot
    % --------------------------------------------------------
    for tt = 1:nT

        target_id = target_des(tt);

        nexttile
        hold on

        if ~cache(tt).has_data
            text(0.5, 0.5, sprintf('Target %d\nnessun dato', target_id), ...
                'Units','normalized', ...
                'HorizontalAlignment','center', ...
                'FontSize',11);
            axis off
            continue;
        end

        t_rel = cache(tt).t_rel;
        psth = cache(tt).psth;

        % cerca il picco inibitorio salvato per questo target/canale
        show_peak = false;
        if has_peak_table && ~isempty(peak_table)
            idx_peak = peak_table.channel_global == ch_global & ...
                       peak_table.target_id == target_id;

            if any(idx_peak)
                peak_time = peak_table.peak_inh_time(find(idx_peak,1,'first'));
                peak_val  = peak_table.peak_inh_abs(find(idx_peak,1,'first'));
                show_peak = ~(isnan(peak_time) || isnan(peak_val));
            end
        end

        % segnale
        plot(t_rel, psth, 'b-', 'LineWidth', 1.6);

        % marker del picco trovato, se presente
        if show_peak
            plot(peak_time, peak_val, 'ro', ...
                'MarkerSize', 7, ...
                'LineWidth', 1.5, ...
                'MarkerFaceColor', 'r');

            text(peak_time, peak_val, sprintf('  %.2f s', peak_time), ...
                'Color', 'r', 'FontSize', 9, 'VerticalAlignment', 'bottom');
        end

        % baseline
        yline(baseline_val, 'k--', 'LineWidth', 1.2);

        % gaze onset
        xline(0, 'k:', 'LineWidth', 1.0);


        % finestra di riferimento opzionale
        if show_window_lines
            xline(win_start_rel, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
            xline(win_end_rel,   '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
        end

        title(sprintf('Target %d', target_id));
        xlabel('Time from gaze onset (s)');
        ylabel('Hz');

        grid on
        box off
        xlim([TMIN_ALL TMAX_ALL]);
        ylim([YMIN_ALL YMAX_ALL]);

        hold off
    end

    if debug_pause
        disp('Premi un tasto per passare al canale successivo...');
        pause;
        if close_after_pause && isvalid(f)
            close(f)
        end
    end
end

fprintf('\nVisualizzazione completata.\n');