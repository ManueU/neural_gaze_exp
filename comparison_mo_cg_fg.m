clear 
% close all
clc
mat_files = { ...
    'motor_BCI02.mat' ... 
    'controlled_BCI02.mat', ...
    'free-gaze_BCI02.mat'
};

n_sets   = 6;
n_trials = 32*ones(1,n_sets);
bin_size = 0.02;

% Stati
NAME_REACH = "Reach";        
offset_reach = 0.02;        
period_reach = 1.0;         

for d = 1:numel(mat_files) 
    disp(mat_files(d)); 
    ds_name = mat_files{d};
    load(ds_name);

    % Movement window
    idx_reach = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == NAME_REACH); 
    start_reach = offset_reach/bin_size;   
    end_reach = start_reach + period_reach/bin_size;

    %Labels
    Y = [];
    for set = 1:n_sets
        Y = [Y; [data(set).Data(1).Resampled.Target_ID]']; 
    end
    
    % Data 
    j = 1; 
    X_pres = cell(sum(n_trials),1);
    X_reach = cell(sum(n_trials),1);
    for set = 1:n_sets
        for trial = 1:n_trials(set)
            tmp_reach = []; 
            for array = 1:2
                tmp_reach = [tmp_reach, data(set).Data(array).Resampled(trial).Task_states{idx_reach,2}(start_reach:end_reach, :)]; 
            end 
            X_reach{j} = mean(tmp_reach./bin_size,1);
            j = j + 1; 
        end   
    end 
    X_reach = cell2mat(X_reach); 
    X_reach = zscore(X_reach); 
    
    
    % SVM 5-fold CV
    k_fold = 5; 
    n_rep = 300; 
    acc_reach = lda_resampled_cv(X_reach, Y, k_fold, n_rep);

    all_acc{d}   = acc_reach;
end 


%% Figure 
classes   = unique(Y);
n_classes  = numel(classes);
chance = 100/n_classes;

labels = {'Motor-only','Controlled-gaze','Free-gaze'};

figure('Color','w');
boxplot(cell2mat(all_acc)*100, ...
    'Colors',[0.3 0.3 0.3], ...
    'Whisker',1.5, 'Symbol','r+', ...
    'Labels',labels);

ylabel('Accuracy (%)');
ylim([0 100]);
yline(chance,'--','Chance','Color',[0.4 0.4 0.4]);



