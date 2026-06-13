clearvars -except mean_baseline_common
close all
clc


% ============================================================
% FILE
% ============================================================
filename_data  = '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat'; %% TO BE CHANGED BASED ON CONDITION
filename_excel = 'modulating.xlsx';

% ============================================================
% PARAMETRI
% ============================================================
sets_plot = [1, 2, 3, 4, 5, 6];
n_sets = numel(sets_plot);

n_arrays   = 2;
n_channels = 96;
bin_size   = 0.02;
w          = 15;

array_names = ["medial", "lateral"];

% finestra rispetto all'onset di gaze
win_start_rel = -1.0;
win_end_rel   =  1.0;

% ============================================================
% CARICAMENTO DATI NEURALI
% ============================================================
S = load(filename_data);
data = S.data;

TS = data(1).Data(1).Interp(1).Task_states;
state_names   = string(TS(:,1));
state_dur_s   = cellfun(@(x) size(x,1) * bin_size, TS(:,2));
state_onset_s = [0; cumsum(state_dur_s(1:end-1))];
get_onset     = @(name) state_onset_s(find(strcmpi(state_names, name), 1, 'first'));

gaze_time = get_onset("gaze");

if isempty(gaze_time) || isnan(gaze_time)
    error('Fase "gaze" non trovata nei Task_states.');
end

% ============================================================
% LETTURA EXCEL
% ============================================================
T3 = readtable(filename_excel, 'Sheet', 'Gaze_on_target');                      %% TO BE CHANGED BASED ON CONDITION
T2 = readtable(filename_excel, 'Sheet', 'Gaze_on_target_details');              %% TO BE CHANGED BASED ON CONDITION

% ---------- FOGLIO 3: canali modulanti ----------
CH = double(T3{:,1});
MOD = double(T3{:,2});
idx = find(MOD == 1);
modulating_channels = CH(idx);

fprintf('Canali modulanti trovati in Foglio3: %d\n', numel(modulating_channels));

% ---------- FOGLIO 2 ----------
% colonna 1 = CH
% colonna 2 = Target
% colonna 3 = flag modulazione (0/1)
CH     = double(T2{:,1});
Target = double(T2{:,2});
Flag   = double(T2{:,3});

% prendi solo righe target con colonna 3 = 1
is_target = ~isnan(Target);
is_mod    = (Flag == 1);
is_mod_ch = ismember(CH, modulating_channels);

sel = is_target & is_mod & is_mod_ch & ~isnan(CH);

pairs = table();
pairs.channel_global = CH(sel);
pairs.target_id      = Target(sel);

pairs = unique(pairs, 'rows');

if isempty(pairs)
    error('Nessuna coppia (canale,target) trovata usando Foglio3 + Foglio2.');
end

fprintf('Coppie canale-target selezionate: %d\n', height(pairs));

target_des = unique(pairs.target_id)';
fprintf('Target presenti: ');
fprintf('%d ', target_des);
fprintf('\n');

% ============================================================
% OUTPUT
% ============================================================
all_rows = table();

% ============================================================
% LOOP SU COPPIE (canale, target)
% ============================================================
for r = 1:height(pairs)

    ch_global = pairs.channel_global(r);
    target_id = pairs.target_id(r);

    array   = ceil(ch_global / n_channels);
    channel = ch_global - (array - 1) * n_channels;

    if array < 1 || array > n_arrays
        fprintf('Salto ch_global=%d: array non valido\n', ch_global);
        continue;
    end

    if channel < 1 || channel > n_channels
        fprintf('Salto ch_global=%d: canale locale non valido\n', ch_global);
        continue;
    end

    if ch_global > numel(mean_baseline_common)
        fprintf('Salto ch_global=%d: baseline non disponibile\n', ch_global);
        continue;
    end

    baseline_val = mean_baseline_common(ch_global);

    if isnan(baseline_val)
        fprintf('Salto ch_global=%d: baseline NaN\n', ch_global);
        continue;
    end

    fprintf('Analizzo target %d | %s | ch %d (glob %d)\n', ...
        target_id, array_names(array), channel, ch_global);

    % --------------------------------------------------------
    % Trial del target nella condizione gaze-only
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
        fprintf('  Nessun trial disponibile\n');
        continue;
    end

    % --------------------------------------------------------
    % PSTH medio smussato
    % --------------------------------------------------------
    psth = smoothdata(mean(firing_rate, 2), 'gaussian', w);
    t = (0:numel(psth)-1) * bin_size;

    this_win_start = max(gaze_time + win_start_rel, t(1));
    this_win_end   = min(gaze_time + win_end_rel,   t(end));

    idx_win = t >= this_win_start & t <= this_win_end;
    if ~any(idx_win)
        fprintf('  Nessun dato nella finestra\n');
        continue;
    end

    t_win_abs = t(idx_win);
    t_win_rel = t_win_abs - gaze_time;
    y_win     = psth(idx_win);

    % ========================================================
    % PICCO INIBITORIO = SEMPLICE MINIMO NELLA FINESTRA [-1 1]
    % ========================================================
    [peak_inh_abs, idx_min] = min(y_win);

    peak_inh_time = t_win_rel(idx_min);          % tempo relativo a gaze onset
    peak_inh_rel  = peak_inh_abs - baseline_val; % valore rispetto a baseline
    peak_inh_amp  = baseline_val - peak_inh_abs; % ampiezza positiva

    % area sotto baseline nella finestra
    below_baseline = max(baseline_val - y_win, 0);
    inh_area_total = trapz(t_win_rel, below_baseline);

    new_row = table( ...
        target_id, ...
        string(array_names(array)), ...
        channel, ...
        ch_global, ...
        baseline_val, ...
        peak_inh_abs, ...
        peak_inh_rel, ...
        peak_inh_amp, ...
        peak_inh_time, ...
        inh_area_total, ...
        'VariableNames', { ...
        'target_id', ...
        'array_name', ...
        'channel_local', ...
        'channel_global', ...
        'baseline', ...
        'peak_inh_abs', ...
        'peak_inh_rel', ...
        'peak_inh_amp', ...
        'peak_inh_time', ...
        'inh_area_total'});

    all_rows = [all_rows; new_row];

    fprintf('  MIN = %.2f Hz | amp = %.2f Hz | time = %.2f s | area = %.2f Hz*s\n', ...
        peak_inh_abs, peak_inh_amp, peak_inh_time, inh_area_total);
end

% ============================================================
% SUMMARY PER TARGET
% ============================================================
fprintf('\n============================================\n');
fprintf('SUMMARY PER TARGET\n');
fprintf('============================================\n');

for i = 1:numel(target_des)
    tid = target_des(i);
    Tt = all_rows(all_rows.target_id == tid, :);

    if isempty(Tt)
        fprintf('Target %d: nessun dato\n', tid);
        continue;
    end

    n_tot = height(Tt);
    med_amp  = median(Tt.peak_inh_amp, 'omitnan');
    med_area = median(Tt.inh_area_total, 'omitnan');

    fprintf('Target %d | coppie analizzate = %d | median amp = %.2f Hz | median area = %.2f Hz*s\n', ...
        tid, n_tot, med_amp, med_area);
end

% ============================================================
% SALVATAGGIO
% ============================================================
save('inhibitory_analysis_got.mat', ...
    'all_rows', 'pairs', 'modulating_channels', ...
    'win_start_rel', 'win_end_rel', 'gaze_time');

disp('File salvato: inhibitory_analysis_got.mat');

% ============================================================
% PREPARAZIONE PLOT
% ============================================================
n_targets = numel(target_des);

n_tot       = zeros(1, n_targets);
median_amp  = nan(1, n_targets);
median_area = nan(1, n_targets);

for i = 1:n_targets
    tid = target_des(i);
    Tt = all_rows(all_rows.target_id == tid, :);

    n_tot(i) = height(Tt);

    if ~isempty(Tt)
        median_amp(i)  = median(Tt.peak_inh_amp, 'omitnan');
        median_area(i) = median(Tt.inh_area_total, 'omitnan');
    end
end

% ============================================================
% FIGURA
% ============================================================
figure('Color','w','Position',[100 100 1400 420]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% ------------------------------------------------------------
% 1) Numero di coppie analizzate per target
% ------------------------------------------------------------
nexttile
bar(target_des, n_tot, 0.7);
xlabel('Target');
ylabel('Count');
title('Responsive channels');
xticks(target_des);
grid on
box off

for i = 1:n_targets
    text(target_des(i), n_tot(i)+0.5, sprintf('n=%d', n_tot(i)), ...
        'HorizontalAlignment','center', 'FontSize',9);
end

% ------------------------------------------------------------
% 2) Boxplot ampiezza del minimo per target
% ------------------------------------------------------------
nexttile
if ~isempty(all_rows)
    boxplot(all_rows.peak_inh_amp, all_rows.target_id, 'Labels', string(target_des));
    hold on
    for i = 1:n_targets
        if ~isnan(median_amp(i))
            plot(i, median_amp(i), 'kd', 'MarkerFaceColor','k', 'MarkerSize',5);
        end
    end
    hold off

    xlabel('Target');
    ylabel('Peak amplitude relative to baseline (Hz)');
    title('Peak inhibitory amplitude');
    grid on
    box off
else
    text(0.5, 0.5, 'No data found', ...
        'Units','normalized', 'HorizontalAlignment','center');
    axis off
end

% ------------------------------------------------------------
% 3) Boxplot area sotto baseline per target
% ------------------------------------------------------------
nexttile
if ~isempty(all_rows)
    boxplot(all_rows.inh_area_total, all_rows.target_id, 'Labels', string(target_des));
    hold on
    for i = 1:n_targets
        if ~isnan(median_area(i))
            plot(i, median_area(i), 'kd', 'MarkerFaceColor','k', 'MarkerSize',5);
        end
    end
    hold off

    xlabel('Target');
    ylabel('Area below baseline (Hz*s)');
    title('Area below baseline');
    grid on
    box off
else
    text(0.5, 0.5, 'No data found', ...
        'Units','normalized', 'HorizontalAlignment','center');
    axis off
end

sgtitle('Gaze-on-target'); %% TO BE CHANGED BASED ON CONDITION