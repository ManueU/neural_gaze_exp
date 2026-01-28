clearvars -except mean_baseline std_baseline
% close all
clc

sets = [2,4,5,6]; 
n_sets = numel(sets); 
n_arrays = 2;
n_trials = 32;
bin_size = 0.02;
period_pre = 1.0;
period_post = 0.5;

filename = { ...
    '../00_Data_extraction/free-gaze_BCI02.mat' ...
    '../00_Data_extraction/motor_BCI02.mat' ...
    '../00_Data_extraction/controlled_BCI02.mat' ...
};

nCond = numel(filename);

%% Costruzione vettore Y e matrice X per SVM
X_all   = cell(nCond,1);
Y_all   = cell(nCond,1);

for d = 1:nCond
    trial_idx = 0;
    fprintf('\nDataset: %s\n', filename{d}); 
    load(filename{d});

    [~, baseName, ext] = fileparts(filename{d});
    ds_name = [baseName ext];

    if strcmp(ds_name, 'controlled_BCI02.mat')
        PRE  = "Gaze";
        POST = "Reach";
    else
        PRE  = "Pres12";
        POST = "Reach";
    end

    Y = [];
    for set_ = 1:n_sets
        set = sets(set_);
        Y = [Y; [data(set).Data(1).Interp.Target_ID]'];
    end
    classes = unique(Y);

    idx_pres  = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == PRE);
    idx_reach = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == POST);

    n_bins_pre  = round(period_pre/bin_size);
    n_bins_post = round(period_post/bin_size);

    start_pre = size(data(1).Data(1).Interp(1).Task_states{idx_pres,2}, 1) - n_bins_pre + 1;
    end_post  = n_bins_post;

    j = 1;
    X = cell(n_trials*n_sets,1);
    for set_ = 1:n_sets
        trial_idx = trial_idx + 1; 
        set = sets(set_);
        for trial = 1:n_trials

            tmp_pre = [];
            tmp_post = [];

            for array = 1:n_arrays
                tmp_pre  = [tmp_pre,  data(set).Data(array).Interp(trial).Task_states{idx_pres,2}(start_pre:end, :)];
                tmp_post = [tmp_post, data(set).Data(array).Interp(trial).Task_states{idx_reach,2}(1:end_post, :)];
            end

            matrix = [tmp_pre; tmp_post];
            firing_rate = matrix./bin_size;
            zscored = (firing_rate - mean_baseline{1,d}(trial_idx, :))./std_baseline{1,d}(trial_idx, :);
            % zscored(isnan(zscored) | isinf(zscored)) = 0;

            X{j} = mean(firing_rate,1);
            j = j + 1;

        end
    end

    X = cell2mat(X);

    % Salvo X e Y per questa condizione
    X_all{d,1} = X;
    Y_all{d,1} = Y;

end

Y_tran = Y_all;         
offset = 0;             
for i = 1:numel(Y_all)
    Y_tran{i} = Y_all{i} + offset;   
    offset = offset + max(Y_all{i});  
end
Y_final = cell2mat(Y_tran);
X_final = cell2mat(X_all);

k_fold = 10;
[acc, cm] = svm_cv(X_final, Y_final, k_fold);


%% Figure (1) - confusion matrix

figure('Color','w');
condShort = {'FG','GC','GT'};  
nTargets  = 8;
classnames = cell(24,1);
for c = 1:nCond      
    for t = 1:nTargets
        idx = (c-1)*nTargets + t;       
        classnames{idx} = sprintf('%s-%d', condShort{c}, t);
    end
end

cc = confusionchart(cm, classnames, ...
    'Normalization','row-normalized', ...
    'RowSummary','off', ...
    'ColumnSummary','off');

cc.GridVisible = 'off';
cc.Title = sprintf('Overall accuracy: %.2f%%', acc*100);
cc.XLabel = '\bf Predicted';
cc.YLabel = '\bf True';

%% Figure (2)
cm_norm = cm ./ sum(cm,2) * 100;  % percentuale

condShort = {'FG','GC','GT'};
nTargets = 8;

classnames = cell(24,1);
for c = 1:3
    for t = 1:nTargets
        idx = (c-1)*nTargets + t;
        classnames{idx} = sprintf('%s-%d', condShort{c}, t);
    end
end

figure('Color','w');
h = heatmap(classnames, classnames, cm_norm);

h.Colormap = parula;
h.ColorLimits = [0 100];
h.CellLabelFormat = '%.1f';  
h.CellLabelColor = 'black'; 
h.XDisplayLabels = classnames;
h.YDisplayLabels = classnames;
h.XLabel = 'Predicted';
h.YLabel = 'True';
h.Title = sprintf('Overall accuracy: %.2f%%', acc*100);

ax = struct(h).Axes;
hold(ax, 'on')

% Linee di separazione
bounds = [8.5 16.5];
for b = bounds
    plot(ax, [b b], [0.5 24.5], 'w-', 'LineWidth', 2);
    plot(ax, [0.5 24.5], [b b], 'w-', 'LineWidth', 2);
end

hold(ax, 'off')


%% Figure (3)
figure('Color','w');
imagesc(cm_norm);
colormap(parula);
colorbar;
axis equal tight;

clear set
set(gca, 'XTick', 1:24, 'XTickLabel', repmat(1:8, 1, 3));
set(gca, 'YTick', 1:24, 'YTickLabel', repmat(1:8, 1, 3));
xlabel('Predicted');
ylabel('True');
title(sprintf('Overall accuracy: %.2f%%', acc*100));

hold on;

% Linee che dividono i blocchi
bounds = [8.5 16.5];
for b = bounds
    plot([b b],[0.5 24.5],'k-','LineWidth',2);
    plot([0.5 24.5],[b b],'k-','LineWidth',2);
end


text(-2, 4.5,  'FG', 'Rotation',90, 'HorizontalAlignment','center', 'FontSize',14, 'FontWeight','bold');
text(-2, 12.5, 'GC', 'Rotation',90, 'HorizontalAlignment','center', 'FontSize',14, 'FontWeight','bold');
text(-2, 20.5, 'GT', 'Rotation',90, 'HorizontalAlignment','center', 'FontSize',14, 'FontWeight','bold');
hold off;
