% =========================================================
% DESCRIZIONE:
%  Script per eseguire una decodifica temporale del Target_ID
%  tramite SVM, utilizzando il firing rate medio estratto da
%  finestre temporali scorrevoli su dati neurali.
%
%  Funzionamento generale:
%    - viene caricato il dataset
%    - si definisce una finestra scorrevole di ampiezza w_length
%      con overlap pari al 50%
%
%  Per ciascuna finestra scorrevole:
%    - per ogni set e per ogni trial vengono estratti i dati
%      da entrambi gli array e concatenati lungo la dimensione
%      dei canali
%    - all’interno della finestra (N_w bin) vengono calcolati
%      i firing rate medi per ciascun canale (Hz)
%    - si costruisce così la matrice delle feature X
%      (tutti i trial × tutti i canali)
%    - si costruisce il vettore delle etichette Y (Target_ID)
%    - si esegue una k-fold cross-validation (k = n) tramite
%      la funzione esterna svm_cv per addestrare un SVM multiclass
%      e stimare l’accuratezza della finestra corrente
%    - vengono salvate:
%         * l’accuratezza complessiva
%         * l’accuratezza per ciascuna classe (diag normalizzata)
%
%  Visualizzazione finale:
%    - andamento temporale dell’accuratezza complessiva e delle
%      accuratezze per classe (curve smussate nel tempo)
%    - indicazione degli eventi/stati del task tramite linee
%      verticali e relative etichette
%    - indicazione del livello di chance (1 / n_classes)
%    - asse temporale basato sul centro di ciascuna finestra
%
%  Parametri principali modificabili:
%    - n_sets, n_arrays, n_trials
%    - bin_size, w_length, overlap
%    - funzione svm_cv (tipo di SVM, k-fold, ecc.)
% =========================================================


clearvars
% close all
clc 

n_sets = 6;
n_arrays = 2;
n_trials = 32;
bin_size = 0.02;

filename = '../00_Data_extraction/gaze_BCI02_withtracker_exclUpdated.mat';
load(filename)

%% Decoding over time with SVM
% bins per trial
N = size(data(1).Data(1).Interp(1).Trial, 1); 
rec_duration = N*bin_size; 

% finestra scorrevole
w_length = 0.6; 
overlap = 0.5*w_length;
N_w = round(w_length/bin_size);
N_o = round(overlap/bin_size);

% etichette Y
Y = []; 
for set = 1:n_sets
    idx = [data(set).Data(1).Interp.Excluded] == 0;
    Y = [Y; [data(set).Data(1).Interp(idx).Target_ID]'];  
end  
classes = unique(Y,'stable');
n_classes = numel(classes);

% numero finestre
n_acc = floor((N - N_w)/(N_w-N_o)) + 1; 
acc_overall = zeros(n_acc,1);
recall_class   = zeros(n_acc, n_classes);
balacc_overall = zeros(n_acc,1);

% loop sulle finestre
start_w = 1; 
end_w = start_w + N_w - 1; 
n_valid = sum( arrayfun(@(s) sum([data(s).Data(1).Interp.Excluded] == 0), 1:n_sets) );

for w = 1:n_acc
    X = cell(n_valid, 1);
    j = 1; 
    for set = 1:n_sets
        for trial = 1:n_trials
            SVM_matrix = []; 
            if data(set).Data(1).Interp(trial).Excluded == 0
                for array = 1:n_arrays
                    SVM_matrix = [SVM_matrix, data(set).Data(array).Interp(trial).Trial(start_w:end_w, :)];   
                end
                X{j} = mean(SVM_matrix./bin_size,1);
                j = j + 1; 
            end        
        end   
    end 
    X = cell2mat(X); 
    
    % k-fold SVM
    % --- Controllo sbilanciamento e k-fold adattivo
    Ycat = categorical(Y);
    counts = countcats(Ycat);
    minCount = min(counts);
    
    % k-fold non può superare il numero di campioni della classe più rara
    k_fold = min(4, minCount);   % oppure min(5, minCount)
    
    if k_fold < 2
        error('Una o più classi hanno <2 campioni: impossibile fare k-fold CV.');
    end
    
    fprintf('Class counts: '); fprintf('%d ', counts); fprintf('\n');
    fprintf('Using k_fold = %d (min class count = %d)\n', k_fold, minCount);

    [acc_overall(w), cm{w}, metrics{w}] = svm_cv(X, Y, k_fold);
    
    balacc_overall(w) = metrics{w}.balancedAccuracy;   % [0..1]

    % accuracy per classe (diag della cm normalizzata per riga) 
    cm_norm = cm{w} ./ max(sum(cm{w},2),1);
    recall_class(w,:) = diag(cm_norm)*100;

    % prossima finestra
    start_w = start_w + (N_w - N_o); 
    end_w = start_w + N_w - 1;  
end 

%% Figure (single panel: overall + per-class recall + correct chance)
figure('Color','White'); hold on;

events_time_tmp = [];
for i = 1:length(data(1).Data(2).Interp(1).Task_states)
    events_time = [events_time_tmp; size(data(1).Data(2).Interp(1).Task_states{i,2},1)*bin_size];
    events_time_tmp = events_time;
end
increment_times = cumsum(events_time);

t = (((0:n_acc-1)*(N_w - N_o)) + N_w/2) * bin_size;

w_smooth = 5;
acc_smooth = smoothdata(acc_overall*100, 'gaussian', w_smooth);
bal_smooth = smoothdata(balacc_overall*100, 'gaussian', w_smooth);

recall_class_smooth = zeros(size(recall_class));
for c = 1:n_classes
    recall_class_smooth(:,c) = smoothdata(recall_class(:,c), 'gaussian', w_smooth);
end

colors = [
    0.839, 0.153, 0.157;
    0.122, 0.467, 0.706;
    0.172, 0.627, 0.172;
    0.580, 0.404, 0.741;
    1.000, 0.498, 0.055;
    0.737, 0.741, 0.133;
    0.549, 0.337, 0.294;
    0.890, 0.466, 0.760;
];

% If more classes than colors, repeat (robust)
if n_classes > size(colors,1)
    colors = repmat(colors, ceil(n_classes/size(colors,1)), 1);
end

for c = 1:n_classes
    plot(t, recall_class_smooth(:,c), ...
        'LineWidth', 0.9, ...
        'Color', colors(c,:), ...
        'HandleVisibility','off');
end

p1 = plot(t, acc_smooth, 'k', 'LineWidth', 1.5, 'DisplayName','Accuracy');
p2 = plot(t, bal_smooth, 'k--', 'LineWidth', 1.5, 'DisplayName','Balanced accuracy');

% Chance levels 
Ycat = categorical(Y);
p = countcats(Ycat) / numel(Ycat);
chance_acc = max(p) * 100;        % majority baseline for accuracy
chance_bal = (1/n_classes) * 100; % uniform chance baseline

% yline(chance_acc, '--', 'Chance (majority)', 'HandleVisibility','off');
yline(chance_bal, ':',  'Chance',      'HandleVisibility','off');

% Events (lines + labels) 
if exist('increment_times','var') && ~isempty(increment_times)
    
    % Free-gaze
    % labels = {"Target cue", "Go cue"}; 
    % xline(increment_times(1:2), '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');

    % Gaze-on-center
    % labels = {"", "Target cue", "Go cue"}; 
    % xline(increment_times(2:3), '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');

    % Gaze-on-target
    labels = {"", "Target cue", "Go cue - gaze", "Go cue - cursor"}; 
    xline(increment_times(2:4), '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');

    % Gaze-only
    % labels = {"", "Target cue", "Go cue"}; 
    % xline(increment_times(2:3), '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');


    ylim([0 100]);
    ax = gca;
    y_pos = ax.YLim(2) - 5;

    for i = 2:4
        x_pos = increment_times(i) - 0.2;
        text(x_pos - 0.05, y_pos, labels{i}, ...
            'HorizontalAlignment','right', ...
            'VerticalAlignment','top', ...
            'Rotation',90, ...
            'FontSize',11);
    end
end

% Axes formatting
ylim([0 100]);
xlim([0 rec_duration]);
yticks(0:20:100);
xlabel('Time (s)');
ylabel('Performance (%)');
box off;

% Legend only for overall curves
legend([p1 p2], 'Location','best');
