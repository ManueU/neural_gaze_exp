% =========================================================
% DESCRIZIONE:
%  Script per eseguire una decodifica temporale del Target_ID
%  tramite SVM utilizzando il firing rate medio estratto da
%  finestre temporali scorrevoli su dati neurali, e per
%  confrontare la overall accuracy tra diverse condizioni
%  sperimentali.
%
%  Funzionamento generale:
%    - per ciascun file .mat elencato in 'filename':
%        * viene caricato il dataset corrispondente
%        * si definisce una finestra scorrevole di durata w_length
%          con overlap pari al 50%
%        * per ogni finestra vengono estratti i dati dai due array,
%          concatenati lungo i canali e convertiti in firing rate medio
%        * viene costruita la matrice delle feature X
%          (tutti i trial × tutti i canali)
%        * viene costruito il vettore Y contenente i Target_ID
%        * si esegue una k-fold cross-validation (k = 5) tramite
%          la funzione esterna svm_cv, ottenendo l’accuratezza
%          per ogni finestra temporale (overall + per classe)
%        * vengono estratti gli stati del task (Task_states) e la
%          loro durata per identificare il tempo di inizio dello
%          stato "Reach"
%        * l’andamento temporale dell’accuratezza viene quindi
%          riallineato rispetto all’onset di Reach
%
%  Output per ciascuna condizione:
%    - overall accuracy nel tempo (smussata)
%    - tempo relativo (t = 0 all’inizio di Reach)
%    - livello di chance (1 / n_classes)
%
%  Visualizzazione finale:
%    - singola figura con le curve di overall accuracy di tutte
%      le condizioni, allineate a t = 0 = onset dello stato reach
%    - indicazione del livello di chance
%
%  Parametri principali modificabili:
%    - n_sets, n_arrays, n_trials
%    - bin_size, w_length, overlap
%    - elenco dei file .mat in 'filename'
%    - nomi delle condizioni in 'cond_names'
%    - funzione svm_cv (kernel, k-fold, ecc.)
% =========================================================


clearvars
close all
clc 

n_sets = 6;
n_arrays = 2;
n_trials = 32;
bin_size = 0.02;

filename = {'../00_Data_extraction/free-gaze_BCI02.mat', ...
            '../00_Data_extraction/motor_BCI02.mat', ...
            '../00_Data_extraction/controlled_BCI02.mat'};

cond_names = { ...
    'Free-gaze', ...
    'Gaze-on-center' ...
    'Gaze-on-target'
    };
 
%% Decoding over time with SVM
acc_class_files = cell(numel(filename), 1);
acc_overall_files = cell(numel(filename), 1);           % non smussata
acc_overall_smooth_files = cell(numel(filename), 1);    % smussata (%)
t_rel_files = cell(numel(filename), 1);                 % tempo relativo a Reach
chance_files = zeros(numel(filename), 1);
for d = 1:numel(filename) 
    disp(filename(d)); 
    load(filename{d});

    % bins per trial
    N = size(data(1).Data(1).Resampled(1).Trial, 1); 
    rec_duration = N*bin_size; 
    
    % finestra scorrevole
    w_length = 0.6; 
    overlap = 0.5*w_length;
    N_w = round(w_length/bin_size);
    N_o = round(overlap/bin_size);
    
    % etichette Y
    Y = []; 
    for set = 1:n_sets
        Y = [Y; [data(set).Data(1).Resampled.Target_ID]']; 
    end  
    classes = unique(Y,'stable');
    n_classes = numel(classes);
    
    % numero finestre
    n_acc = floor((N - N_w)/(N_w-N_o)) + 1; 
    acc_overall = zeros(n_acc,1);
    acc_class = zeros(n_acc, n_classes);
    cm = cell(n_acc, 1);
    
    % loop sulle finestre
    start_w = 1; 
    end_w = start_w + N_w - 1; 
    
    for w = 1:n_acc
        X = cell(n_trials*n_sets,1);
        j = 1; 
        for set = 1:n_sets
            for trial = 1:n_trials
                SVM_matrix = []; 
                for array = 1:n_arrays
                    SVM_matrix = [SVM_matrix, data(set).Data(array).Resampled(trial).Trial(start_w:end_w, :)]; 
                end 
                X{j} = mean(SVM_matrix./bin_size,1);
                j = j + 1; 
            end   
        end 
        X = cell2mat(X); 
        
        % k-fold SVM
        k_fold = 5; 
        [acc_overall(w), cm{w}] = svm_cv(X, Y, k_fold);
        
        % accuracy per classe (diag della cm normalizzata per riga) 
        cm_norm = cm{w} ./ max(sum(cm{w},2),1);
        acc_class(w,:) = diag(cm_norm)*100;
    
        % prossima finestra
        start_w = start_w + (N_w - N_o); 
        end_w = start_w + N_w - 1;  
    end 

    % Ricavo tempi degli stati ed onset di Reach
    Task_states = data(1).Data(2).Resampled(1).Task_states;
    n_states = size(Task_states, 1);
    events_time = zeros(n_states, 1);
    for i = 1:n_states
        events_time(i) = size(Task_states{i, 2}, 1) * bin_size;
    end
    increment_times = cumsum(events_time);       
    labels = string(Task_states(:, 1)); 
    start_states = [0; increment_times(1:end-1)];
    idxReach = find(labels == "Reach", 1, 'first');
    t0 = start_states(idxReach);   
    
    t = (((0:n_acc-1) * (N_w - N_o)) + N_w/2) * bin_size;  % tempo assoluto (centro finestra)
    t_rel = t - t0;                                        % tempo relativo a onset Reach

    w_smooth = 5;
    acc_smooth_overall = smoothdata(acc_overall, 'gaussian', w_smooth) * 100; % in %

    acc_class_files{d} = acc_class;
    acc_overall_files{d} = acc_overall;
    acc_overall_smooth_files{d} = acc_smooth_overall;
    t_rel_files{d} = t_rel;
    chance_files(d) = (1 / n_classes) * 100;
    
end 

%% Figure
figure('Color', 'white'); hold on;

n_cond = numel(filename);
cols = [
    0.5  0.5  0.5;   % grigio
    0.855, 0.537, 0.600;   % rosa
    0.2  0.85 0.75;  % verde acqua
];

for d = 1:n_cond
    plot(t_rel_files{d}, acc_overall_smooth_files{d}, ...
        'LineWidth', 1.5, ...
        'Color', cols(d, :), ...
        'DisplayName', cond_names{d});
end

% livello di chance
chance_level = mean(chance_files);
yline(chance_level, '--k', 'Chance', 'HandleVisibility', 'off');

% linea verticale per l'onset di Reach (t = 0)
xline(0, ':k', 'Reach onset', 'LabelVerticalAlignment', 'top', ...
      'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');

xlabel('Time from reach onset (s)');
ylabel('Overall accuracy (%)');
legend('Location', 'best');

% limiti in x basati su tutte le condizioni
all_t = cell2mat(t_rel_files(:)');
xlim([min(all_t), max(all_t)]);
ylim([0, 100]);