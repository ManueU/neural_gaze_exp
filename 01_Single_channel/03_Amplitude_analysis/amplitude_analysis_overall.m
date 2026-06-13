clearvars
close all
clc

% ============================================================
% LOAD
% ============================================================
load('amplitude_analysis_global_baseline.mat', 'all_rows');
load('responsive_channels_free_gaze.mat', 'responsive_matrix');

if ~exist('all_rows','var') || isempty(all_rows)
    error('La tabella all_rows non è disponibile o è vuota.');
end

if ~exist('responsive_matrix','var') || isempty(responsive_matrix)
    error('responsive_matrix non disponibile o vuota.');
end

% ============================================================
% PARAMETRI
% ============================================================
n_arrays   = 2;
n_channels = 96;
n_total_channels = n_arrays * n_channels;

array_names = ["medial", "lateral"];

% ============================================================
% CANALI RESPONSIVI IN FREE-GAZE (regardless target)
% ============================================================
responsive_free_gaze_any = any(responsive_matrix, 1);
free_gaze_channels = find(responsive_free_gaze_any);
n_free = numel(free_gaze_channels);

fprintf('\n============================================\n');
fprintf('CANALI RESPONSIVI IN FREE-GAZE\n');
fprintf('============================================\n');
fprintf('Canali responsivi in free-gaze regardless target: %d\n', n_free);
fprintf('Percentuale sui %d canali totali: %.1f%%\n', ...
    n_total_channels, 100 * n_free / n_total_channels);

% ============================================================
% TABELLA DI PARTENZA
% ============================================================
T = all_rows;

% Overall combinata: max tra eccitazione e inibizione
T.resp_center_total = max([T.exc_center, T.inh_center_amp], [], 2);
T.resp_target_total = max([T.exc_target, T.inh_target_amp], [], 2);

T.center_present_any = T.exc_center_present | T.inh_center_present;
T.target_present_any = T.exc_target_present | T.inh_target_present;

% ============================================================
% COSTRUZIONE SUMMARY A LIVELLO DI CANALE
% ============================================================
channel_summary_all = table();
channel_summary_exc = table();
channel_summary_inh = table();

for ch_global = 1:n_total_channels

    Tc = T(T.channel_global == ch_global, :);

    array_id = ceil(ch_global / n_channels);
    channel_local = ch_global - (array_id - 1) * n_channels;

    if isempty(Tc)
        row_all = make_empty_row(ch_global, array_id, channel_local, array_names(array_id));
        row_exc = make_empty_row(ch_global, array_id, channel_local, array_names(array_id));
        row_inh = make_empty_row(ch_global, array_id, channel_local, array_names(array_id));

    else
        % ====================================================
        % OVERALL COMBINATA
        % ====================================================
        responds_center_any = any(Tc.center_present_any);
        responds_target_any = any(Tc.target_present_any);
        responds_any = responds_center_any || responds_target_any;

        resp_center_overall = max(Tc.resp_center_total);
        resp_target_overall = max(Tc.resp_target_total);

        stronger_in_center = responds_any && (resp_center_overall > resp_target_overall);
        stronger_in_target = responds_any && (resp_target_overall > resp_center_overall);
        equal_response     = responds_any && (resp_target_overall == resp_center_overall);

        selective_center_only = responds_center_any && ~responds_target_any;
        selective_target_only = responds_target_any && ~responds_center_any;

        row_all = table( ...
            ch_global, array_id, channel_local, string(array_names(array_id)), ...
            responds_any, responds_center_any, responds_target_any, ...
            resp_center_overall, resp_target_overall, ...
            stronger_in_center, stronger_in_target, equal_response, ...
            selective_center_only, selective_target_only, ...
            'VariableNames', { ...
            'channel_global', 'array_id', 'channel_local', 'array_name', ...
            'responds_any', 'responds_center_any', 'responds_target_any', ...
            'resp_center_overall', 'resp_target_overall', ...
            'stronger_in_center', 'stronger_in_target', 'equal_response', ...
            'selective_center_only', 'selective_target_only'});

        % ====================================================
        % SOLO ECCITATORIA
        % ====================================================
        exc_center_present = any(Tc.exc_center_present);
        exc_target_present = any(Tc.exc_target_present);
        exc_any = exc_center_present || exc_target_present;

        exc_center_overall = max(Tc.exc_center);
        exc_target_overall = max(Tc.exc_target);

        exc_stronger_in_center = exc_any && (exc_center_overall > exc_target_overall);
        exc_stronger_in_target = exc_any && (exc_target_overall > exc_center_overall);
        exc_equal_response     = exc_any && (exc_target_overall == exc_center_overall);

        exc_selective_center_only = exc_center_present && ~exc_target_present;
        exc_selective_target_only = exc_target_present && ~exc_center_present;

        row_exc = table( ...
            ch_global, array_id, channel_local, string(array_names(array_id)), ...
            exc_any, exc_center_present, exc_target_present, ...
            exc_center_overall, exc_target_overall, ...
            exc_stronger_in_center, exc_stronger_in_target, exc_equal_response, ...
            exc_selective_center_only, exc_selective_target_only, ...
            'VariableNames', { ...
            'channel_global', 'array_id', 'channel_local', 'array_name', ...
            'responds_any', 'responds_center_any', 'responds_target_any', ...
            'resp_center_overall', 'resp_target_overall', ...
            'stronger_in_center', 'stronger_in_target', 'equal_response', ...
            'selective_center_only', 'selective_target_only'});

        % ====================================================
        % SOLO INIBITORIA
        % ====================================================
        inh_center_present = any(Tc.inh_center_present);
        inh_target_present = any(Tc.inh_target_present);
        inh_any = inh_center_present || inh_target_present;

        inh_center_overall = max(Tc.inh_center_amp);
        inh_target_overall = max(Tc.inh_target_amp);

        inh_stronger_in_center = inh_any && (inh_center_overall > inh_target_overall);
        inh_stronger_in_target = inh_any && (inh_target_overall > inh_center_overall);
        inh_equal_response     = inh_any && (inh_target_overall == inh_center_overall);

        inh_selective_center_only = inh_center_present && ~inh_target_present;
        inh_selective_target_only = inh_target_present && ~inh_center_present;

        row_inh = table( ...
            ch_global, array_id, channel_local, string(array_names(array_id)), ...
            inh_any, inh_center_present, inh_target_present, ...
            inh_center_overall, inh_target_overall, ...
            inh_stronger_in_center, inh_stronger_in_target, inh_equal_response, ...
            inh_selective_center_only, inh_selective_target_only, ...
            'VariableNames', { ...
            'channel_global', 'array_id', 'channel_local', 'array_name', ...
            'responds_any', 'responds_center_any', 'responds_target_any', ...
            'resp_center_overall', 'resp_target_overall', ...
            'stronger_in_center', 'stronger_in_target', 'equal_response', ...
            'selective_center_only', 'selective_target_only'});
    end

    channel_summary_all = [channel_summary_all; row_all];
    channel_summary_exc = [channel_summary_exc; row_exc];
    channel_summary_inh = [channel_summary_inh; row_inh];
end

% ============================================================
% STAMPA RISULTATI
% ============================================================
fprintf('\n============================================\n');
fprintf('OVERALL COMBINATA\n');
fprintf('============================================\n');
print_summary(channel_summary_all, n_total_channels, free_gaze_channels);

fprintf('\n============================================\n');
fprintf('RISPOSTA ECCITATORIA\n');
fprintf('============================================\n');
print_summary(channel_summary_exc, n_total_channels, free_gaze_channels);

fprintf('\n============================================\n');
fprintf('RISPOSTA INIBITORIA\n');
fprintf('============================================\n');
print_summary(channel_summary_inh, n_total_channels, free_gaze_channels);

% ============================================================
% SALVATAGGIO DATI
% ============================================================
% save('overall_channel_characterization_separate_exc_inh.mat', ...
%     'channel_summary_all', 'channel_summary_exc', 'channel_summary_inh', ...
%     'free_gaze_channels', 'responsive_free_gaze_any');

disp('File salvato: overall_channel_characterization_separate_exc_inh.mat');

% ============================================================
% FIGURA COMPLETA: OVERALL + EXC + INH
% Tutte le percentuali sono normalizzate sui canali free-gaze
% ============================================================
figure('Color','w','Position',[100 50 1450 950]);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

% ---------------- OVERALL ----------------
plot_summary_panel(channel_summary_all, free_gaze_channels, ...
    'Responsive channels', 1);
plot_dominance_panel(channel_summary_all, free_gaze_channels, ...
    'Larger response', 2);
plot_selectivity_panel(channel_summary_all, free_gaze_channels, ...
    'Condition selectivity', 3);

% ---------------- ECCITATORIA ----------------
plot_summary_panel(channel_summary_exc, free_gaze_channels, ...
    '', 4);
plot_dominance_panel(channel_summary_exc, free_gaze_channels, ...
    '', 5);
plot_selectivity_panel(channel_summary_exc, free_gaze_channels, ...
    '', 6);

% ---------------- INIBITORIA -----------------
plot_summary_panel(channel_summary_inh, free_gaze_channels, ...
    '', 7);
plot_dominance_panel(channel_summary_inh, free_gaze_channels, ...
    '', 8);
plot_selectivity_panel(channel_summary_inh, free_gaze_channels, ...
    '', 9);

sgtitle('Channel characterization');

saveas(gcf, 'overall_channel_characterization_overall_exc_inh.png');

% ============================================================
% FUNZIONI LOCALI
% ============================================================
function row = make_empty_row(ch_global, array_id, channel_local, array_name)
    row = table( ...
        ch_global, array_id, channel_local, string(array_name), ...
        false, false, false, ...
        0, 0, ...
        false, false, false, ...
        false, false, ...
        'VariableNames', { ...
        'channel_global', 'array_id', 'channel_local', 'array_name', ...
        'responds_any', 'responds_center_any', 'responds_target_any', ...
        'resp_center_overall', 'resp_target_overall', ...
        'stronger_in_center', 'stronger_in_target', 'equal_response', ...
        'selective_center_only', 'selective_target_only'});
end

function print_summary(S, n_total_channels, free_gaze_channels)
    n_resp_any_total   = sum(S.responds_any);
    n_no_resp_total    = n_total_channels - n_resp_any_total;

    Sfg = S(free_gaze_channels, :);
    n_free = numel(free_gaze_channels);

    n_resp_any   = sum(Sfg.responds_any);
    n_no_resp    = n_free - n_resp_any;

    n_str_center = sum(Sfg.stronger_in_center);
    n_str_target = sum(Sfg.stronger_in_target);
    n_equal      = sum(Sfg.equal_response);

    n_sel_center = sum(Sfg.selective_center_only);
    n_sel_target = sum(Sfg.selective_target_only);
    n_both       = sum(Sfg.responds_center_any & Sfg.responds_target_any);

    fprintf('Canali totali: %d\n', n_total_channels);
    fprintf('Canali responsivi in free-gaze: %d\n', n_free);

    fprintf('\nSu tutti i canali:\n');
    fprintf('- canali responsivi: %d\n', n_resp_any_total);
    fprintf('- canali non responsivi: %d\n', n_no_resp_total);

    fprintf('\nSul sottoinsieme free-gaze:\n');
    fprintf('- canali responsivi: %d\n', n_resp_any);
    fprintf('- canali non responsivi: %d\n', n_no_resp);
    fprintf('- percentuale responsivi: %.1f%%\n', 100*n_resp_any/n_free);

    fprintf('\nDominanza di condizione (sui canali responsivi del sottoinsieme free-gaze):\n');
    fprintf('- maggiore in gaze-on-center: %d\n', n_str_center);
    fprintf('- maggiore in gaze-on-target: %d\n', n_str_target);
    fprintf('- uguale: %d\n', n_equal);

    fprintf('\nSelettivita'' di condizione (sui canali responsivi del sottoinsieme free-gaze):\n');
    fprintf('- solo gaze-on-center: %d\n', n_sel_center);
    fprintf('- solo gaze-on-target: %d\n', n_sel_target);
    fprintf('- entrambe le condizioni: %d\n', n_both);

    if n_resp_any > 0
        fprintf('\nPercentuali sui canali responsivi del sottoinsieme free-gaze:\n');
        fprintf('- maggiore in center: %.1f%%\n', 100*n_str_center/n_resp_any);
        fprintf('- maggiore in target: %.1f%%\n', 100*n_str_target/n_resp_any);
        fprintf('- solo center: %.1f%%\n', 100*n_sel_center/n_resp_any);
        fprintf('- solo target: %.1f%%\n', 100*n_sel_target/n_resp_any);
        fprintf('- entrambe: %.1f%%\n', 100*n_both/n_resp_any);
    end
end

function plot_summary_panel(S, free_gaze_channels, ttl, tile_idx)
    nexttile(tile_idx)

    Sfg = S(free_gaze_channels, :);

    n_resp_any = sum(Sfg.responds_any);
    n_no_resp  = numel(free_gaze_channels) - n_resp_any;

    vals = 100 * [n_resp_any, n_no_resp] / numel(free_gaze_channels);
    counts = [n_resp_any, n_no_resp];

    bar(vals, 0.6);
    ylabel('Percentage of channels (%)');
    xticks(1:2);
    xticklabels({'Responding','Non-responding'});
    ylim([0 100]);
    title(ttl);
    grid on
    box off

    for i = 1:numel(vals)
        text(i, vals(i)+4, sprintf('%.1f%% (n=%d)', vals(i), counts(i)), ...
            'HorizontalAlignment','center', 'FontSize',10);
    end
end

function plot_dominance_panel(S, free_gaze_channels, ttl, tile_idx)
    nexttile(tile_idx)

    Sfg = S(free_gaze_channels, :);

    n_resp_any   = sum(Sfg.responds_any);
    n_str_center = sum(Sfg.stronger_in_center);
    n_str_target = sum(Sfg.stronger_in_target);
    n_equal      = sum(Sfg.equal_response);

    counts = [n_str_center, n_str_target, n_equal];

    if n_resp_any > 0
        vals = 100 * counts / n_resp_any;
    else
        vals = [0 0 0];
    end

    bar(vals, 0.6);
    ylabel('Percentage of channels (%)');
    xticks(1:3);
    xticklabels({'Center','Target','Equal'});
    ylim([0 100]);
    title(ttl);
    grid on
    box off

    for i = 1:numel(vals)
        text(i, vals(i)+4, sprintf('%.1f%% (n=%d)', vals(i), counts(i)), ...
            'HorizontalAlignment','center', 'FontSize',10);
    end
end

function plot_selectivity_panel(S, free_gaze_channels, ttl, tile_idx)
    nexttile(tile_idx)

    Sfg = S(free_gaze_channels, :);

    n_resp_any   = sum(Sfg.responds_any);
    n_sel_center = sum(Sfg.selective_center_only);
    n_sel_target = sum(Sfg.selective_target_only);
    n_both       = sum(Sfg.responds_center_any & Sfg.responds_target_any);

    counts = [n_sel_center, n_sel_target, n_both];

    if n_resp_any > 0
        vals = 100 * counts / n_resp_any;
    else
        vals = [0 0 0];
    end

    bar(vals, 0.6);
    ylabel('Percentage of channels (%)');
    xticks(1:3);
    xticklabels({'Center','Target','Both'});
    ylim([0 100]);
    title(ttl);
    grid on
    box off

    for i = 1:numel(vals)
        text(i, vals(i)+4, sprintf('%.1f%% (n=%d)', vals(i), counts(i)), ...
            'HorizontalAlignment','center', 'FontSize',10);
    end
end