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
offset_reach = 0.2;         % inizio a +200 ms dall'onset del Reach
period_reach = 0.2;          % durata 200 ms (-> 200–400 ms)

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

    X_pres = zscore(X_pres); 
    X_reach = zscore(X_reach); 
    
    
    % LDA 5-fold CV
    k_fold = 5; 
    n_rep = 300; 
    acc_pre   = lda_resampled_cv(X_pres,  Y, k_fold, n_rep);
    acc_reach = lda_resampled_cv(X_reach, Y, k_fold, n_rep);

    all_acc{d}   = [acc_pre(:),  acc_reach(:)];
end 


%% Statistical test
acc1 = all_acc{1};
acc2 = all_acc{2};
acc = [acc1(:,1), acc2(:,1), acc1(:,2), acc2(:,2)];

delta_pre = acc(:,1) - acc(:,2); 
delta_reach = acc(:,3) - acc(:,4); 
mean_delta_pre = mean(delta_pre);
mean_delta_reach = mean(delta_reach);

if mean_delta_pre > 0
    p_pre = 2 * (sum(delta_pre < 0) / n_rep);
else
    p_pre = 2 * (sum(delta_pre > 0) / n_rep);
end
fprintf('Mean Δ = %.4f   p_pre = %.4g\n', mean_delta_pre, p_pre);

if mean_delta_reach > 0
    p_reach = 2 * (sum(delta_reach < 0) / n_rep);
else
    p_reach = 2 * (sum(delta_reach > 0) / n_rep);
end
fprintf('Mean Δ = %.4f   p_reach = %.4g\n', mean_delta_reach, p_reach);

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

%% Figure 
classes   = unique(Y);
n_classes  = numel(classes);
chance = 100/n_classes;

grp = [ ...
    repmat({'Prep'},  size(acc1,1),1); ...
    repmat({'Prep'},  size(acc2,1),1); ...
    repmat({'Move'}, size(acc1,1),1); ...
    repmat({'Move'}, size(acc2,1),1)];
ds  = [ ...
    repmat({'Motor-only'}, size(acc1,1),1); ...
    repmat({'Free-gaze'}, size(acc2,1),1); ...
    repmat({'Motor-only'}, size(acc1,1),1); ...
    repmat({'Free-gaze'}, size(acc2,1),1)];

labels = {'Motor-only','Free-gaze','Motor-only','Free-gaze'};

figure('Color','w');
boxplot(acc*100, {grp, ds}, ...
    'Colors',[0.3 0.3 0.3], ...
    'Whisker',1.5, 'Symbol','r+', ...
    'Labels',labels, ...              
    'FactorGap',[12 1], ...
    'FactorSeparator',1);

ylabel('Accuracy (%)');
ylim([0 100]);
yline(chance,'--','Chance','Color',[0.4 0.4 0.4]);


ax = gca;
xt = ax.XTick;             
yl = ax.YLim;              

ytop = yl(2) - 0.08*diff(yl);   
text(mean(xt(1:2)), ytop, 'Preparation', 'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', 'FontWeight','bold');
text(mean(xt(3:4)), ytop, 'Movement', 'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', 'FontWeight','bold');

hold on
y_gap = 4; % spazio verticale tra annotazioni

% Prep
if p_pre < 0.05
    y_pos = yl(2) - 5*y_gap;
    plot([1 2], [y_pos y_pos], 'k', 'LineWidth',1.2);
    text(1.5, y_pos + 1.5, getStars(p_pre), 'HorizontalAlignment','center', 'FontSize',12);
end

% Reach
if p_reach < 0.05
    y_pos = yl(2) - 5*y_gap;
    plot([3.4 4.4], [y_pos y_pos], 'k', 'LineWidth',1.2);
    text(3.9, y_pos + 1.5, getStars(p_reach), 'HorizontalAlignment','center', 'FontSize',12);
end
hold off


