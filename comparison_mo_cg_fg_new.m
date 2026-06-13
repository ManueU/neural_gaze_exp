clear all 
close all 
clc
mat_files = { ...
    'motor_BCI02.mat' ... 
    % 'controlled_BCI02.mat', ...
    'free-gaze_BCI02.mat'
};

n_sets = 6; 
n_trials = 32*ones(1,n_sets);  
bin_size = 0.02;

figure('Color', 'White')
labels_cond = ["Motor-only", "Free-gaze"]; 
colors_cond = [233, 114, 55; 
          59, 125, 35]/255; 
for d = 1:numel(mat_files) 
    disp(mat_files(d)); 
    ds_name = mat_files{d};
    load(ds_name);

    N            = size(data(1).Data(1).Resampled(1).Trial, 1); 
    rec_duration = N*bin_size; 
    w_length     = 0.6; 
    overlap      = 0.5*w_length;
    N_w          = round(w_length/bin_size);
    N_o          = round(overlap/bin_size);
    
    % Labels
    tmp = []; 
    for set = 1:n_sets
        Y = [tmp; [data(set).Data(1).Resampled.Target_ID]']; 
        tmp = Y;
    end 
    classes = unique(Y,'stable');
    n_classes = numel(classes);
    
    n_acc = floor((N - N_w)/(N_w - N_o)) + 1; 
    start_w = 1; 
    end_w = start_w + N_w - 1; 
    
    acc_overall = zeros(n_acc,1);
    acc_class   = zeros(n_acc, n_classes);
    w_middle = zeros(n_acc,1);

    % Loop sulle finestre
    for w = 1:n_acc
        j = 1; 
        X = cell(sum(n_trials),1);
        for set = 1:n_sets
            for trial = 1:n_trials(set)
                SVM_matrix = []; 
                for array = 1:2
                    SVM_matrix = [SVM_matrix, data(set).Data(array).Resampled(trial).Trial(start_w:end_w, :)]; 
                end 
                X{j} = mean(SVM_matrix./bin_size,1);
                j = j + 1; 
            end   
        end 
        X = cell2mat(X); 
        
        k_fold = 5; 
        [acc_overall(w)] = lda_cv(X, Y, k_fold);
        % cm_norm = cm{w} ./ max(sum(cm{w},2),1);
        % acc_class(w,:) = diag(cm_norm)*100;
        w_middle(w) = start_w + N_w/2; 
    
        start_w = start_w + (N_w - N_o); 
        end_w = start_w + N_w - 1; 
    end 

    % Figure     
    events_time_tmp = []; 
    for i = 1:length(data(1).Data(2).Resampled(1).Task_states)
        events_time = [events_time_tmp; size(data(1).Data(2).Resampled(1).Task_states{i,2},1)*bin_size];
        events_time_tmp = events_time; 
    end 
    increment_times = cumsum(events_time); 
    labels = string(data(1).Data(2).Resampled(1).Task_states(:,1));

    idx_reach = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == "Reach"); 
    start_reach = increment_times(idx_reach - 1); 
    
    % figure('Color', 'White')
    colors = [
        0.839, 0.153, 0.157;  % rosso
        0.122, 0.467, 0.706;  % blu
        0.172, 0.627, 0.172;  % verde
        0.580, 0.404, 0.741;  % viola
        1.000, 0.498, 0.055;  % arancione
        0.737, 0.741, 0.133;  % giallo oliva
        0.549, 0.337, 0.294;  % marrone
        0.890, 0.466, 0.760;  % rosa
    ];
    alpha = 0.5;
    
    t = (((0:n_acc-1)*(N_w - N_o)) + N_w/2)*bin_size - start_reach;
    idx_time = find(t >= -0.8 & t <= 0.8);
    % t = (((0:n_acc-1)*(N_w - N_o)) + N_w/2)*bin_size;
    % for c = 1:n_classes
    %     acc_smooth = smoothdata(acc_class(:, c), 'gaussian', 4);
    %     plot(t, acc_smooth, 'LineWidth', 1.0, 'Color', [colors(c,:),  alpha], 'HandleVisibility','off'), hold on
    % end
    acc_smooth_overall = smoothdata(acc_overall, 'gaussian', 4);
    % plot(t(idx_time), acc_smooth_overall(idx_time)*100, 'LineWidth', 1.5, 'Color', 'k', 'DisplayName', 'Overall'), hold on
    plot(t(idx_time), acc_smooth_overall(idx_time)*100, 'LineWidth', 1.5, 'Color', colors_cond(d,:), 'DisplayName',labels_cond(d)), hold on


    % if exist('increment_times','var') && ~isempty(increment_times)
    %    xline(increment_times, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
    %    k = min(numel(labels), numel(increment_times));
    %    x_times = [0; increment_times(:)];
    %    x_text  = x_times(1:k) + diff(x_times(1:k+1))/2;
         % ylim([0 100]);
    %    ax = gca; 
    %    y_text = (ax.YLim(2) - 10)*ones(1,length(x_text)); 
         % text(x_text, y_text, labels(1:k), 'HorizontalAlignment','center');
    % end
    % yline((1/n_classes)*100,'-', 'Chance', 'HandleVisibility','off'); 
    % legend show; 
    % xlabel('Time (s)');
    % ylabel('Accuracy (%)');
    % xlim([0 rec_duration]);

    acc_overall_files{d} = [w_middle, acc_overall]; 
    % acc_class_files{d} = acc_class;  
end 
time_vect = -0.6:0.2:0.6; 
xticks(time_vect);
xline(0, '--', 'HandleVisibility','off');
yline((1/n_classes)*100,'-', 'Chance', 'HandleVisibility','off'); 
xlim([-0.6 0.6]);
ylim([0 100]);
legend show; 
xlabel('Time (s)');
ylabel('Overall accuracy (%)');



%% Comparison   
clc
clearvars -except acc_overall_files acc_class_files mat_files bin_size n_trials n_sets

for d = 1:numel(mat_files) 
    disp(mat_files(d)); 
    ds_name = mat_files{d};
    load(ds_name);

    % Ricerca del massimo per definire la finestra
    delta_t = 0.1; 
    id = find(acc_overall_files{1,d}(:,2) == max(acc_overall_files{1,d}(:,2))); 
    N_peak = acc_overall_files{1,d}(id,1); 
    start_reach = N_peak - round(delta_t/bin_size);   
    end_reach = N_peak + round(delta_t/bin_size);

    %Labels
    Y = [];
    for set = 1:n_sets
        Y = [Y; [data(set).Data(1).Resampled.Target_ID]']; 
    end
    
    % Data 
    j = 1; 
    X_reach = cell(sum(n_trials),1);
    for set = 1:n_sets
        for trial = 1:n_trials(set)
            tmp_reach = []; 
            for array = 1:2
                tmp_reach = [tmp_reach, data(set).Data(array).Resampled(trial).Trial(start_reach:end_reach, :)]; 
            end 
            X_reach{j} = mean(tmp_reach./bin_size,1);
            j = j + 1; 
        end   
    end 
    X_reach = cell2mat(X_reach); 
    % X_reach = zscore(X_reach); 
    
    
    % SVM 5-fold CV
    k_fold = 5; 
    n_rep = 300; 
    acc_reach = svm_resampled_cv(X_reach, Y, k_fold, n_rep);

    all_acc{d}   = acc_reach;
end 

%% Statistical analysis
acc1 = all_acc{1};
acc2 = all_acc{2};
acc3 = all_acc{3};

delta_12 = acc2 - acc1; 
delta_13 = acc3 - acc1; 
delta_23 = acc3 - acc2;
mean_delta_12 = mean(delta_12);
mean_delta_13 = mean(delta_13);
mean_delta_23 = mean(delta_23);

if mean_delta_12 > 0
    p_12 = 2 * (sum(delta_12 < 0) / n_rep);
else
    p_12 = 2 * (sum(delta_12 > 0) / n_rep);
end
fprintf('Mean Δ = %.4f   p_12 = %.4g\n', mean_delta_12, p_12);

if mean_delta_13 > 0
    p_13 = 2 * (sum(delta_13 < 0) / n_rep);
else
    p_13 = 2 * (sum(delta_13 > 0) / n_rep);
end
fprintf('Mean Δ = %.4f   p_13 = %.4g\n', mean_delta_13, p_13);

if mean_delta_23 > 0
    p_23 = 2 * (sum(delta_23 < 0) / n_rep);
else
    p_23 = 2 * (sum(delta_23 > 0) / n_rep);
end
fprintf('Mean Δ = %.4f   p_23 = %.4g\n', mean_delta_23, p_23);


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

labels = {'Motor-only','Controlled-gaze','Free-gaze'};

figure('Color','w');
boxplot(cell2mat(all_acc)*100, ...
    'Colors',[0.3 0.3 0.3], ...
    'Whisker',1.5, 'Symbol','r+', ...
    'Labels',labels);

ylabel('Accuracy (%)');
ylim([0 100]);
yline(chance,'--','Chance','Color',[0.4 0.4 0.4]);


