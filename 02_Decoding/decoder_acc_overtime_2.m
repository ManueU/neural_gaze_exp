clear all 
close all 
clc
load('gaze_BCI03_2.mat')

n_sets = 6; 
n_trials = 32*ones(1,n_sets);  
% n_trials = [40,32,32,32];  
bin_size = 0.02;

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
    [acc_overall(w), cm{w}] = svm_cv(X, Y, k_fold);
    cm_norm = cm{w} ./ max(sum(cm{w},2),1);
    acc_class(w,:) = diag(cm_norm)*100;

    start_w = start_w + (N_w - N_o); 
    end_w = start_w + N_w - 1; 
end 

%% Figure
events_time_tmp = []; 
for i = 1:length(data(1).Data(2).Resampled(1).Task_states)
    events_time = [events_time_tmp; size(data(1).Data(2).Resampled(1).Task_states{i,2},1)*bin_size];
    events_time_tmp = events_time; 
end 
increment_times = cumsum(events_time); 
labels = string(data(1).Data(2).Resampled(1).Task_states(:,1));

figure('Color', 'White')
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

t = (((0:n_acc-1)*(N_w - N_o)) + N_w/2)*bin_size;
for c = 1:n_classes
    acc_smooth = smoothdata(acc_class(:, c), 'gaussian', 4);
    plot(t, acc_smooth, 'LineWidth', 1.0, 'Color', [colors(c,:),  alpha], 'HandleVisibility','off'), hold on
end
acc_smooth_overall = smoothdata(acc_overall, 'gaussian', 4);
plot(t, acc_smooth_overall*100, 'LineWidth', 1.5, 'Color', 'k', 'DisplayName','Overall'), hold on

if exist('increment_times','var') && ~isempty(increment_times)
   xline(increment_times, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
   k = min(numel(labels), numel(increment_times));
   x_times = [0; increment_times(:)];
   x_text  = x_times(1:k) + diff(x_times(1:k+1))/2;
   ylim([0 100]);
   ax = gca; 
   y_text = (ax.YLim(2) - 10)*ones(1,length(x_text)); 
   text(x_text, y_text, labels(1:k), 'HorizontalAlignment','center');
end
yline((1/n_classes)*100,'-', 'Chance', 'HandleVisibility','off'); 
legend show; 
xlabel('Time (s)');
ylabel('Accuracy (%)');
xlim([0 rec_duration]);

