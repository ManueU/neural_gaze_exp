clear 
% close all
clc
mat_files = { ...
    'motor_BCI02.mat' ... 
    'free-gaze_BCI02.mat', ...         
};

n_sets   = 6;
n_trials = 32*ones(1,n_sets);
bin_size = 0.02;

% Stati
NAME_PRE   = "Pres12";       % finestra preparatoria (ultimi 200 ms dello stato Pres12)
NAME_REACH = "Reach";        % finestra movimento (200–400 ms dall'inizio dello stato Reach)
period_pre   = 0.2;          % ultimi 200 ms del Pres12
offset_reach = 0.2;          % inizio a +200 ms dall'onset del Reach
period_reach = 0.4;          % durata 200 ms (-> 200–400 ms)

for d = 1:numel(mat_files) 
    ds_name = mat_files{d};
    load(ds_name);

    % Preparatory window
    idx_pres = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == NAME_PRE); 
    start_pres = size(data(1).Data(2).Resampled(1).Task_states{idx_pres,2}, 1) - (period_pre/bin_size); 
    end_pres = size(data(1).Data(2).Resampled(1).Task_states{idx_pres,2}, 1); 
    
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
            tmp_pres = []; 
            tmp_reach = []; 
            for array = 1:2
                tmp_pres = [tmp_pres, data(set).Data(array).Resampled(trial).Task_states{idx_pres,2}(start_pres:end_pres, :)]; 
                tmp_reach = [tmp_reach, data(set).Data(array).Resampled(trial).Task_states{idx_reach,2}(start_reach:end_reach, :)]; 
            end 
            X_pres{j} = mean(tmp_pres./bin_size,1);
            X_reach{j} = mean(tmp_reach./bin_size,1);
            j = j + 1; 
        end   
    end 
    X_pres = cell2mat(X_pres); 
    X_reach = cell2mat(X_reach); 
    
    
    % LDA 5-fold CV
    k_fold = 5; 
    acc_pre   = lda_cv(X_pres,   Y, k_fold);
    acc_reach = lda_cv(X_reach, Y, k_fold);

    all_acc{d}   = [acc_pre(:),  acc_reach(:)];
end 

%% Figure 
classes   = unique(Y);
n_classes  = numel(classes);
chance = 100/n_classes;

acc1 = all_acc{1};
acc2 = all_acc{2};

acc = [acc1(:,1); acc2(:,1); acc1(:,2); acc2(:,2)];
grp = [ ...
    repmat({'Prep'},  size(acc1,1),1); ...
    repmat({'Prep'},  size(acc2,1),1); ...
    repmat({'Reach'}, size(acc1,1),1); ...
    repmat({'Reach'}, size(acc2,1),1)];
ds  = [ ...
    repmat({'Motor-only'}, size(acc1,1),1); ...
    repmat({'Free-gaze'}, size(acc2,1),1); ...
    repmat({'Motor-only'}, size(acc1,1),1); ...
    repmat({'Free-gaze'}, size(acc2,1),1)];

figure('Color','w');
boxplot(acc*100, {grp, ds}, ...
    'colors', [0.3 0.3 0.3], ...
    'whisker', 1.5, 'symbol', 'r+');

ylabel('Accuracy (%)');
ylim([0 100]);
yline(chance,'--','Chance','Color',[0.4 0.4 0.4]); 

%% Statistical test
acc_prep_motor  = acc1(:,1); % Motor-only, Prep
acc_prep_free   = acc2(:,1); % Free-gaze, Prep
acc_reach_motor = acc1(:,2); % Motor-only, Reach
acc_reach_free  = acc2(:,2); % Free-gaze, Reach

% Test di normalità (Shapiro-Wilk o Lilliefors)
[h1_prep_motor, p1_prep_motor]   = lillietest(acc_prep_motor);
[h1_prep_free, p1_prep_free]     = lillietest(acc_prep_free);
[h1_reach_motor, p1_reach_motor] = lillietest(acc_reach_motor);
[h1_reach_free, p1_reach_free]   = lillietest(acc_reach_free);

if ~any([h1_prep_motor h1_prep_free h1_reach_motor h1_reach_free])
    % t-test a due campioni
    [~, p_prep]  = ttest2(acc_prep_motor, acc_prep_free);
    [~, p_reach] = ttest2(acc_reach_motor, acc_reach_free);
    test_used = 't-test (independent samples)';
else
    % Test non parametrico
    p_prep  = signrank(acc_prep_motor, acc_prep_free);
    p_reach = signrank(acc_reach_motor, acc_reach_free);
    test_used = 'Wilcoxon rank-sum';
end

function stars = getStars(p)
    if p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    elseif p < 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end
end
hold on
yl = ylim;
y_gap = 4; % spazio verticale tra annotazioni

% Prep
if p_prep < 0.05
    y_pos = yl(2) - y_gap;
    plot([1 2], [y_pos y_pos], 'k', 'LineWidth',1.2);
    text(1.5, y_pos + 1.5, getStars(p_prep), 'HorizontalAlignment','center', 'FontSize',12);
end

% Reach
if p_reach < 0.05
    y_pos = yl(2) - 2*y_gap;
    plot([3 4], [y_pos y_pos], 'k', 'LineWidth',1.2);
    text(3.5, y_pos + 1.5, getStars(p_reach), 'HorizontalAlignment','center', 'FontSize',12);
end
hold off
fprintf('\nStatistical comparison (%s)\n', test_used);
fprintf('\nPrep:  p = %.4f (%s)\n', p_prep, getStars(p_prep));
fprintf('Reach: p = %.4f (%s)\n', p_reach, getStars(p_reach));
