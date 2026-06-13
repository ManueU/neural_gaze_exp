clearvars -except mean_baseline std_baseline
% close all
clc

sets_corr = [2,4,5,6];
n_sets = numel(sets_corr); 
n_arrays = 2;
n_targets = 8; 
n_trials = 32; 
n_channels = 96;
bin_size = 0.02;
period_pre = 0.1; 
period_post = 0.5; 

filename = { 
    '../00_Data_extraction/free-gaze_BCI02.mat'...
    '../00_Data_extraction/motor_BCI02.mat',...
    '../00_Data_extraction/controlled_BCI02.mat'
    };

%% Data construction
condition_matrix = [];
for d = 1:numel(filename)
    load(filename{d}); 

    [~, baseName, ext] = fileparts(filename{d});
    ds_name = [baseName ext];

    if strcmp(ds_name, 'controlled_BCI02.mat')
        PRE = "Gaze";
        POST = "Reach";
    elseif strcmp(ds_name, 'gaze_BCI02.mat')
        PRE = "Pres12";
        POST = "Gaze";
    else
        PRE = "Pres12";
        POST = "Reach";
    end

    idx_pres = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == PRE); 
    idx_reach = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == POST); 
    if isempty(idx_pres) || isempty(idx_reach)
        error('Stati PRE/POST non trovati: controlla PRE/POST');
    end

    pre_bins  = max(1, round(period_pre/bin_size));
    post_bins = max(1, round(period_post/bin_size));
    correlation_matrix = [];
    for array = 1:n_arrays
        mean_across_time = nan(n_targets, n_channels);
        for channel = 1:n_channels
            ch_global = (array-1)*n_channels + channel;
            for target = 1:n_targets
                firing_rate = []; 
                for set_ = 1:n_sets
                    set = sets_corr(set_);
                    idx = find([data(set).Data(array).Interp.Target_ID] == target); 
                   
                    for j = 1:length(idx) 
                        tmp_pre = data(set).Data(array).Interp(idx(j)).Task_states{idx_pres, 2}(end-pre_bins+1:end, channel); 
                        tmp_post = data(set).Data(array).Interp(idx(j)).Task_states{idx_reach,2}(1:post_bins, channel); 
                        vect = [tmp_pre; tmp_post];  
                        firing_rate = [firing_rate, vect./bin_size];  
                    end
                end
                zscored = (mean(firing_rate,2) - mean_baseline.by_targets{1,d}(target, ch_global))./std_baseline.by_targets{1,d}(target, ch_global);
                mean_across_time(target,channel) = mean(zscored); 
            end 
        end 
        correlation_matrix = [correlation_matrix, mean_across_time]; 
    end 
    condition_matrix = [condition_matrix; correlation_matrix];
    condition(d) = {correlation_matrix};
end 

%% Correlation 
R = corr(condition_matrix', 'Rows','pairwise');


%% Figure (1)
gaze_names = {'FG','GC','GT'};  
n_gaze = numel(gaze_names);
C = n_targets*n_gaze;

assert(all(size(R) == [C C]), 'R must be CxC with C=n_targets*n_gaze');

figure('Color','w'); 
imagesc(R);
axis image;
colormap(parula);
colorbar;
hold on;

for g = 1:n_gaze-1
    x = g*n_targets + 0.5;
    line([x x], [0.5 C+0.5], 'Color','k', 'LineWidth', 1.8);
    line([0.5 C+0.5], [x x], 'Color','k', 'LineWidth', 1.8);
end

xticks(1:C);
yticks(1:C);

xlab = strings(1,C);
ylab = strings(1,C);
idx = 0; 
for d = 1:n_gaze
    for t = 1:n_targets
        idx = idx + 1; 
        xlab(idx) = t;
        ylab(idx) = t;
    end
end
xticklabels(xlab);
yticklabels(ylab);

for g = 1:n_gaze
    center = (g-1)*n_targets + (n_targets+1)/2;
    text(center, 0.0, gaze_names{g}, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontWeight','bold', 'Interpreter','none', 'Clipping','off');
    text(-1.5, center, gaze_names{g}, ...
        'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
        'FontWeight','bold', 'Interpreter','none', 'Clipping','off');
end

clear set
set(gca, 'TickLength',[0 0], 'FontSize',10);
hold off;
