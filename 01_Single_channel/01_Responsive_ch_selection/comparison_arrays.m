%% Responsive channel overlap map

% Visualizes the spatial distribution of responsive channels across the
% three experimental conditions and highlights their overlap on the array.
% Each electrode is classified according to the number of conditions in
% which it exhibits significant modulation.

clc
clear

array = 1;

% Load channel map
load('ChannelMap_BCI02.mat')
motor_medial  = ChannelMap.ChannelNumbers{1,1};
motor_lateral = ChannelMap.ChannelNumbers{1,3};
motor_electrodes = {motor_medial, motor_lateral};

motor_map = motor_electrodes{array};

% Load responsive channels for all 3 experimental conditions
load('responsive_channels_free_gaze.mat')
ch_cond1 = unique(results.channel_global);

load('responsive_channels_gaze_on_center.mat')
ch_cond2 = unique(results.channel_global);

load('responsive_channels_gaze_on_target.mat')
ch_cond3 = unique(results.channel_global);

% Mask of valid positions
valid = ~isnan(motor_map);

% Logic matrices
mod1 = false(size(motor_map));
mod2 = false(size(motor_map));
mod3 = false(size(motor_map));

mod1(valid) = ismember(motor_map(valid), ch_cond1);
mod2(valid) = ismember(motor_map(valid), ch_cond2);
mod3(valid) = ismember(motor_map(valid), ch_cond3);

% Number of conditions in which each channel is responsive
overlap = double(mod1) + double(mod2) + double(mod3);

% Convert to double and assign NaN to missing electrode positions for plotting
mod1_plot = double(mod1);
mod2_plot = double(mod2);
mod3_plot = double(mod3);
overlap_plot = overlap;

mod1_plot(~valid) = NaN;
mod2_plot(~valid) = NaN;
mod3_plot(~valid) = NaN;
overlap_plot(~valid) = NaN;

% Figure
figure('Color','w','Position',[100 100 1000 250])

% --- Condition 1 ---
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
        
        % up
        if i == 1 || ~valid(i-1,j)
            plot([j-0.5 j+0.5],[i-0.5 i-0.5],'k','LineWidth',1)
        end
        
        % down
        if i == nr || ~valid(i+1,j)
            plot([j-0.5 j+0.5],[i+0.5 i+0.5],'k','LineWidth',1)
        end
        
        % left
        if j == 1 || ~valid(i,j-1)
            plot([j-0.5 j-0.5],[i-0.5 i+0.5],'k','LineWidth',1)
        end
        
        % right
        if j == nc || ~valid(i,j+1)
            plot([j+0.5 j+0.5],[i-0.5 i+0.5],'k','LineWidth',1)
        end
        
    end
end

hold off

% --- Condition 2 ---
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
        
        % up
        if i == 1 || ~valid(i-1,j)
            plot([j-0.5 j+0.5],[i-0.5 i-0.5],'k','LineWidth',1)
        end
        
        % down
        if i == nr || ~valid(i+1,j)
            plot([j-0.5 j+0.5],[i+0.5 i+0.5],'k','LineWidth',1)
        end
        
        % left
        if j == 1 || ~valid(i,j-1)
            plot([j-0.5 j-0.5],[i-0.5 i+0.5],'k','LineWidth',1)
        end
        
        % right
        if j == nc || ~valid(i,j+1)
            plot([j+0.5 j+0.5],[i-0.5 i+0.5],'k','LineWidth',1)
        end
        
    end
end

hold off

% --- Condition 3 ---
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
        
        % up
        if i == 1 || ~valid(i-1,j)
            plot([j-0.5 j+0.5],[i-0.5 i-0.5],'k','LineWidth',1)
        end
        
        % down
        if i == nr || ~valid(i+1,j)
            plot([j-0.5 j+0.5],[i+0.5 i+0.5],'k','LineWidth',1)
        end
        
        % left
        if j == 1 || ~valid(i,j-1)
            plot([j-0.5 j-0.5],[i-0.5 i+0.5],'k','LineWidth',1)
        end
        
        % right
        if j == nc || ~valid(i,j+1)
            plot([j+0.5 j+0.5],[i-0.5 i+0.5],'k','LineWidth',1)
        end
        
    end
end

hold off

% --- Overlap ---
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
        
        % up
        if i == 1 || ~valid(i-1,j)
            plot([j-0.5 j+0.5],[i-0.5 i-0.5],'k','LineWidth',1)
        end
        
        % down
        if i == nr || ~valid(i+1,j)
            plot([j-0.5 j+0.5],[i+0.5 i+0.5],'k','LineWidth',1)
        end
        
        % left
        if j == 1 || ~valid(i,j-1)
            plot([j-0.5 j-0.5],[i-0.5 i+0.5],'k','LineWidth',1)
        end
        
        % right
        if j == nc || ~valid(i,j+1)
            plot([j+0.5 j+0.5],[i-0.5 i+0.5],'k','LineWidth',1)
        end
        
    end
end

hold off

[nr,nc] = size(overlap_plot);

% Add numbers only to gray cells (overlap = 1 or 2)
for i = 1:nr
    for j = 1:nc
        
        if isnan(overlap_plot(i,j))
            continue
        end

        % Labels
        labels = [];
        if mod1(i,j), labels = [labels, 1]; end
        if mod2(i,j), labels = [labels, 2]; end
        if mod3(i,j), labels = [labels, 3]; end

        if ~isempty(labels) && (overlap(i,j) == 1 || overlap(i,j) == 2)

            txt = strjoin(string(labels), ',');

            % text color
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