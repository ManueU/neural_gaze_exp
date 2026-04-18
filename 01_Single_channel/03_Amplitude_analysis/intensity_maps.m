clc
clear
close all

%% =========================================================
% PARAMETRI
% ==========================================================
array_to_plot = 2;   % 1 = medial, 2 = lateral

cond_names = {'Free-gaze', 'Gaze-on-center', 'Gaze-on-target'};
result_files = { ...
    'responsive_channels_free_gaze.mat', ...
    'responsive_channels_gaze_on_center.mat', ...
    'responsive_channels_gaze_on_target.mat'};

use_percentile_scaling = true;   % se true usa percentile per evitare outlier
perc_val = 95;

%% =========================================================
% CARICAMENTO MAPPA CANALI
% ==========================================================
load('ChannelMap_BCI02.mat')

motor_medial  = ChannelMap.ChannelNumbers{1,1};
motor_lateral = ChannelMap.ChannelNumbers{1,3};
motor_electrodes = {motor_medial, motor_lateral};

motor_map = motor_electrodes{array_to_plot};
valid = ~isnan(motor_map);

[nr, nc] = size(motor_map);

%% =========================================================
% PREALLOCAZIONE
% 3 condizioni x 2 tipi risposta
% ==========================================================
exc_maps = cell(1,3);
inh_maps = cell(1,3);

all_exc_vals = [];
all_inh_vals = [];

%% =========================================================
% LOOP SULLE CONDIZIONI
% ==========================================================
for c = 1:3

    load(result_files{c}, 'results');

    % numero canali globali massimo trovato nel file
    max_ch = max(results.channel_global);

    % vettori con intensità massima per canale
    exc_strength = nan(max_ch,1);
    inh_strength = nan(max_ch,1);

    % canali presenti nel file
    ch_list = unique(results.channel_global);

    for ch = ch_list'

        idx = results.channel_global == ch;

        % massimo picco eccitatorio relativo alla baseline
        if any(idx)
            exc_vals = results.peak_exc_rel(idx);
            inh_vals = results.peak_inh_rel(idx);

            if ~all(isnan(exc_vals))
                exc_strength(ch) = max(exc_vals, [], 'omitnan');
            end

            if ~all(isnan(inh_vals))
                inh_strength(ch) = max(abs(inh_vals), [], 'omitnan');
            end
        end
    end

    % proiezione sulla mappa fisica dell'array
    exc_map = nan(size(motor_map));
    inh_map = nan(size(motor_map));

    for i = 1:nr
        for j = 1:nc
            if ~valid(i,j)
                continue
            end

            ch = motor_map(i,j);

            if ch <= numel(exc_strength)
                exc_map(i,j) = exc_strength(ch);
            end

            if ch <= numel(inh_strength)
                inh_map(i,j) = inh_strength(ch);
            end
        end
    end

    exc_maps{c} = exc_map;
    inh_maps{c} = inh_map;

    all_exc_vals = [all_exc_vals; exc_map(valid)];
    all_inh_vals = [all_inh_vals; inh_map(valid)];
end

%% =========================================================
% SCALE COLORI SEPARATE
% confrontabili tra condizioni, ma separate per tipo risposta
% ==========================================================

% ----- eccitatoria
exc_vals_all = all_exc_vals(~isnan(all_exc_vals));

if isempty(exc_vals_all)
    clim_exc = [0 1];
else
    if use_percentile_scaling
        vmax_exc = prctile(exc_vals_all, perc_val);
    else
        vmax_exc = max(exc_vals_all);
    end
    if vmax_exc == 0
        vmax_exc = 1;
    end
    clim_exc = [0 vmax_exc];
end

% ----- inibitoria
inh_vals_all = all_inh_vals(~isnan(all_inh_vals));

if isempty(inh_vals_all)
    clim_inh = [0 1];
else
    if use_percentile_scaling
        vmax_inh = prctile(inh_vals_all, perc_val);
    else
        vmax_inh = max(inh_vals_all);
    end
    if vmax_inh == 0
        vmax_inh = 1;
    end
    clim_inh = [0 vmax_inh];
end

%% =========================================================
% DIFFERENZA: gaze-on-center vs gaze-on-target
% c=2 --> gaze-on-center
% c=3 --> gaze-on-target
% ==========================================================
diff_exc_abs = abs(exc_maps{3} - exc_maps{2});
diff_inh_abs = abs(inh_maps{3} - inh_maps{2});

all_diff_vals = [diff_exc_abs(valid); diff_inh_abs(valid)];
all_diff_vals = all_diff_vals(~isnan(all_diff_vals));

if isempty(all_diff_vals)
    diff_clim = [0 1];
else
    if use_percentile_scaling
        diff_max = prctile(all_diff_vals, perc_val);
    else
        diff_max = max(all_diff_vals);
    end
    if diff_max == 0
        diff_max = 1;
    end
    diff_clim = [0 diff_max];
end

%% =========================================================
% COLORMAP RISPOSTE: bianco sporco -> arancio -> rosso
% ==========================================================
n = 256;
x = linspace(0,1,n)';

cmap_resp = [ ...
    0.98 + 0.02*x, ...
    0.98 - 0.65*x, ...
    0.96 - 0.96*x.^1.4 ...
];

%% =========================================================
% FIGURA 1: ECCITATORIA
% ==========================================================
figure('Color','w','Position',[100 100 1300 700])

for c = 1:3
    subplot(1,3,c)
    imagesc(exc_maps{c}, 'AlphaData', ~isnan(exc_maps{c}))
    axis image off
    title(cond_names{c})
    colormap(gca, cmap_resp)
    caxis(clim_exc)
    colorbar
    hold on
    draw_array_outline(valid)
    hold off
end

sgtitle('Excitatory response', 'FontWeight', 'bold')

%% =========================================================
% FIGURA 2: INIBITORIA
% ==========================================================
figure('Color','w','Position',[100 100 1300 700])

for c = 1:3
    subplot(1,3,c)
    imagesc(inh_maps{c}, 'AlphaData', ~isnan(inh_maps{c}))
    axis image off
    title(cond_names{c})
    colormap(gca, cmap_resp)
    caxis(clim_inh)
    colorbar
    hold on
    draw_array_outline(valid)
    hold off
end

sgtitle('Inhibitory response', 'FontWeight', 'bold')

%% =========================================================
% FIGURA 3: DIFFERENZA ASSOLUTA CENTER vs TARGET
% stessa colormap delle figure principali
% ==========================================================
figure('Color','w','Position',[100 100 1000 500])

n = 256;
x = linspace(0,1,n)';

cmap_diff = [ ...
    0.98 - 0.60*x, ...   % R ↓
    0.98 - 0.40*x, ...   % G ↓ (più lento → azzurro)
    0.96 - 0.10*x ...    % B quasi alto → blu soft
];

subplot(1,2,1)
imagesc(diff_exc_abs, 'AlphaData', ~isnan(diff_exc_abs))
axis image off
title('Gaze-on-target - Gaze-on-center')
colormap(gca, cmap_diff)
caxis(diff_clim)
colorbar
hold on
draw_array_outline(valid)
hold off

subplot(1,2,2)
imagesc(diff_inh_abs, 'AlphaData', ~isnan(diff_inh_abs))
axis image off
title('Gaze-on-target - Gaze-on-center')
colormap(gca, cmap_diff)
caxis(diff_clim)
colorbar
hold on
draw_array_outline(valid)
hold off


%% =========================================================
% FUNZIONE LOCALE: contorno array
% ==========================================================
function draw_array_outline(valid)

[nr,nc] = size(valid);

for i = 1:nr
    for j = 1:nc

        if ~valid(i,j)
            continue
        end

        % sopra
        if i == 1 || ~valid(i-1,j)
            plot([j-0.5 j+0.5],[i-0.5 i-0.5],'k','LineWidth',1)
        end

        % sotto
        if i == nr || ~valid(i+1,j)
            plot([j-0.5 j+0.5],[i+0.5 i+0.5],'k','LineWidth',1)
        end

        % sinistra
        if j == 1 || ~valid(i,j-1)
            plot([j-0.5 j-0.5],[i-0.5 i+0.5],'k','LineWidth',1)
        end

        % destra
        if j == nc || ~valid(i,j+1)
            plot([j+0.5 j+0.5],[i-0.5 i+0.5],'k','LineWidth',1)
        end
    end
end

end