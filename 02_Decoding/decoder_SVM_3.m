clear 
close all
clc

mat_files = { ...
    'free-gaze_BCI02.mat' ... 
    'motor_BCI02.mat'...
    'controlled_BCI02.mat'...
    'gaze_BCI02.mat'
};

n_sets   = 6; 
n_trials = 32*ones(1,n_sets);  
bin_size = 0.02;

for d = 1:numel(mat_files) 
    disp(mat_files(d)); 
    ds_name = mat_files{d};
    load(ds_name);

    % Time windows and states names
    if strcmp(ds_name, 'gaze_BCI02.mat')
        NAME_PRE     = "Pres12"; 
        NAME_REACH   = "Gaze"; 
        period_pre   = 0.1; 
        period_reach = 0.1; 
    else 
        if strcmp(ds_name, 'controlled_BCI02.mat')
            NAME_PRE     = "Gaze"; 
            NAME_REACH   = "Reach"; 
            period_pre   = 0.1; 
            period_reach = 0.5;
        else 
            NAME_PRE     = "Pres12"; 
            NAME_REACH   = "Reach"; 
            period_pre   = 0.1; 
            period_reach = 0.5;
        end 
    end 
    
    % Preparatory window
    idx_pres = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == NAME_PRE); 
    start_pres = size(data(1).Data(2).Resampled(1).Task_states{idx_pres,2}, 1) - (period_pre/bin_size); 
    
    % Movement window
    idx_reach = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == NAME_REACH); 
    end_reach = period_reach/bin_size;
    
    % Labels
    tmp = []; 
    for set = 1:n_sets
        Y = [tmp; [data(set).Data(1).Resampled.Target_ID]']; 
        tmp = Y;
    end 
    
    % Data
    j = 1; 
    X = cell(sum(n_trials),1);
    for set = 1:n_sets
        for trial = 1:n_trials(set)
            tmp_pres = []; 
            tmp_reach = []; 
            for array = 1:2
                tmp_pres = [tmp_pres, data(set).Data(array).Resampled(trial).Task_states{idx_pres,2}(start_pres:end, :)]; 
                tmp_reach = [tmp_reach, data(set).Data(array).Resampled(trial).Task_states{idx_reach,2}(1:end_reach, :)]; 
            end 
            matrix = [tmp_pres; tmp_reach]; 
            X{j} = mean(matrix./bin_size,1);
            j = j + 1; 
        end   
    end 
    X = cell2mat(X); 
    
    k_fold = 5; 
    [acc, cm] = svm_cv(X, Y, k_fold); 
    figure('Color','w');
    classes = unique(Y);
    classnames = arrayfun(@(c) sprintf('Target %d', c), classes, 'UniformOutput', false);
    cc = confusionchart(cm, classnames, 'Normalization','row-normalized', 'RowSummary','off', 'ColumnSummary','off');
    cc.Title  = sprintf('Confusion Matrix - Location Decoder (Accuracy: %.2f%%)', acc*100);
    cc.XLabel = 'Predicted Target';
    cc.YLabel = 'True Target';

    cm_all{d,1} = cm;
    acc_all{d,1} = acc; 
end


%% Barplot
figure('Color','w');
acc = cell2mat(acc_all)*100;   % accuracy in percentuale
nClassi = numel(unique(Y));    % numero di target/classi (per la chance level)
chance = (1/nClassi)*100;      % livello di chance in %

b = bar(acc, 0.5, 'FaceColor', [0.2 0.4 0.7], 'EdgeColor', 'none'); 
ylim([0 100]);
ylabel('Performance (%)');
xlabel('Condition');
xticks(1:4);
xticklabels({'Free-gaze', 'Motor-only', 'Controlled-gaze', 'Gaze-only'});
title('SVM overall accuracy across conditions');

% --- Linea di chance level ---
hold on;
yline(chance, '-k', sprintf('Chance level (%.1f%%)', chance), ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', ...
    'FontSize', 9, 'LineWidth', 1.2);
grid on;
box off;

