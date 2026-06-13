clc
clear

array = 1;

% Carica mappa canali
load('ChannelMap_BCI02.mat')
motor_medial  = ChannelMap.ChannelNumbers{1,1};
motor_lateral = ChannelMap.ChannelNumbers{1,3};
motor_electrodes = {motor_medial, motor_lateral};

motor_map = motor_electrodes{array};

% Carica canali responsive nelle 3 condizioni
load('responsive_channels_free_gaze.mat')
ch_cond1 = unique(results.channel_global);

load('responsive_channels_gaze_on_center.mat')
ch_cond2 = unique(results.channel_global);

load('responsive_channels_gaze_on_target.mat')
ch_cond3 = unique(results.channel_global);

% Maschera posizioni valide
valid = ~isnan(motor_map);

% Matrici logiche
mod1 = false(size(motor_map));
mod2 = false(size(motor_map));
mod3 = false(size(motor_map));

mod1(valid) = ismember(motor_map(valid), ch_cond1);
mod2(valid) = ismember(motor_map(valid), ch_cond2);
mod3(valid) = ismember(motor_map(valid), ch_cond3);

% Numero di condizioni in cui ogni canale è modulante
overlap = double(mod1) + double(mod2) + double(mod3);

% Per il plotting: converti a double e metti NaN dove non esistono elettrodi
mod1_plot = double(mod1);
mod2_plot = double(mod2);
mod3_plot = double(mod3);
overlap_plot = overlap;

mod1_plot(~valid) = NaN;
mod2_plot(~valid) = NaN;
mod3_plot(~valid) = NaN;
overlap_plot(~valid) = NaN;

% Figura
figure('Color','w','Position',[100 100 1000 250])

%% --- CONDITION 1 ---
subplot(1,4,1)
imagesc(mod1_plot)
axis image off
title('Free-gaze')
colormap(gca, [1 1 1; 0 0 0])
caxis([0 1])

hold on

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

hold off

%% --- CONDITION 2 ---
subplot(1,4,2)
imagesc(mod2_plot)
axis image off
title('Gaze-on-center')
colormap(gca, [1 1 1; 0 0 0])
caxis([0 1])

hold on

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

hold off

%% --- CONDITION 3 ---
subplot(1,4,3)
imagesc(mod3_plot)
axis image off
title('Gaze-on-target')
colormap(gca, [1 1 1; 0 0 0])
caxis([0 1])

hold on

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

hold off

%% --- OVERLAP ---
subplot(1,4,4)
imagesc(overlap_plot)
axis image off
title('Overlap')
colormap(gca, [1 1 1; 0.75 0.75 0.75; 0.35 0.35 0.35; 0 0 0])
caxis([0 3])

hold on

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

hold off

[nr,nc] = size(overlap_plot);

% aggiungi numeri SOLO nelle celle grigie (overlap = 1 o 2)
for i = 1:nr
    for j = 1:nc
        
        if isnan(overlap_plot(i,j))
            continue
        end

        % costruisci label condizioni
        labels = [];
        if mod1(i,j), labels = [labels, 1]; end
        if mod2(i,j), labels = [labels, 2]; end
        if mod3(i,j), labels = [labels, 3]; end

        if ~isempty(labels) && (overlap(i,j) == 1 || overlap(i,j) == 2)

            txt = strjoin(string(labels), ',');

            % colore testo
            txt_color = 'w';

            text(j, i, txt, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontSize',8, ...
                'FontWeight','bold', ...
                'Color', txt_color);
        end
    end
end

%% --- LEGEND ---
hold on

% colori (uguali alla colormap overlap)
c_none = [1 1 1];
c_1    = [0.75 0.75 0.75];
c_2    = [0.35 0.35 0.35];
c_3    = [0 0 0];

h1 = patch(nan, nan, c_none);
h2 = patch(nan, nan, c_1);
h3 = patch(nan, nan, c_2);
h4 = patch(nan, nan, c_3);

legend([h1 h2 h3 h4], ...
    {'No modulation','Modulation in 1 cond','Modulation in 2 cond','Modulation'}, ...
    'Location','best', ...
    'Orientation','vertical', ...
    'Box','off');