clc
clear
% close all

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
                inh_strength(ch) = min(inh_vals, [], 'omitnan');
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
    clim_inh = [-1 0];
else
    if use_percentile_scaling
        vmin_inh = prctile(inh_vals_all, 100 - perc_val);
    else
        vmin_inh = min(inh_vals_all);
    end
    if vmin_inh == 0
        vmin_inh = -1;
    end
    clim_inh = [vmin_inh 0];
end

%% =========================================================
% DIFFERENZA SIGNED NORMALIZZATA: gaze-on-target vs gaze-on-center
% misura = (target - center) / max(target, center)
% c=2 --> gaze-on-center
% c=3 --> gaze-on-target
% ==========================================================

center_exc = exc_maps{2};
target_exc = exc_maps{3};

center_inh = inh_maps{2};
target_inh = inh_maps{3};

max_exc = max(target_exc, center_exc);
max_inh = max(target_inh, center_inh);

diff_exc_norm = (target_exc - center_exc) ./ max_exc;
diff_inh_norm = (target_inh - center_inh) ./ max_inh;

% evita Inf dove il denominatore è zero
diff_exc_norm(~isfinite(diff_exc_norm)) = NaN;
diff_inh_norm(~isfinite(diff_inh_norm)) = NaN;

% scala colori simmetrica attorno a zero
all_diff_vals = [diff_exc_norm(valid); diff_inh_norm(valid)];
all_diff_vals = all_diff_vals(~isnan(all_diff_vals));

if isempty(all_diff_vals)
    diff_clim = [-1 1];
else
    if use_percentile_scaling
        diff_lim = prctile(abs(all_diff_vals), perc_val);
    else
        diff_lim = max(abs(all_diff_vals));
    end
    if diff_lim == 0
        diff_lim = 1;
    end
    diff_clim = [-diff_lim diff_lim];
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

sgtitle('Excitatory response (Hz)', 'FontWeight', 'bold')

%% =========================================================
% FIGURA 2: INIBITORIA
% ==========================================================
figure('Color','w','Position',[100 100 1300 700])

for c = 1:3
    subplot(1,3,c)
    imagesc(inh_maps{c}, 'AlphaData', ~isnan(inh_maps{c}))
    axis image off
    title(cond_names{c})
    colormap(gca, flipud(cmap_resp))
    caxis(clim_inh)
    colorbar
    hold on
    draw_array_outline(valid)
    hold off
end

sgtitle('Inhibitory response (Hz)', 'FontWeight', 'bold')

%% =========================================================
% FIGURA 3: DIFFERENZA SIGNED NORMALIZZATA
% (target - center) / media(target, center)
% colormap divergente centrata su zero
% ==========================================================
figure('Color','w','Position',[100 100 1000 500])

n = 256;
half = n/2;

% ---- blu -> bianco
x1 = linspace(0,1,half)';
neg = [ ...
    0.10 + 0.90*x1, ...   % R
    0.25 + 0.75*x1, ...   % G
    0.75 + 0.25*x1 ...    % B
];

% ---- bianco -> ocra
x2 = linspace(0,1,half)';
pos = [ ...
    1.00 - 0.20*x2, ...   % R (leggermente meno di bianco)
    1.00 - 0.45*x2, ...   % G
    1.00 - 0.80*x2 ...    % B (crolla → giallo/arancio)
];

cmap_diff = [neg; pos];

subplot(1,2,1)
imagesc(diff_exc_norm, 'AlphaData', ~isnan(diff_exc_norm))
axis image off
title('Excitatory response')
colormap(gca, cmap_diff)
caxis(diff_clim)
cb1 = colorbar;
hold on
draw_array_outline(valid)
hold off

subplot(1,2,2)
imagesc(diff_inh_norm, 'AlphaData', ~isnan(diff_inh_norm))
axis image off
title('Inhibitory response')
colormap(gca, cmap_diff)
caxis(diff_clim)
cb2 = colorbar;
hold on
draw_array_outline(valid)
hold off

sgtitle({ ...
    'Gaze-on-target vs Gaze-on-center', ...
    '(signed difference normalized to max)'}, ...
    'FontWeight', 'bold')

