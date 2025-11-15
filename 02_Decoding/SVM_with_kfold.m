% =========================================================
% DESCRIZIONE:
%  Script per addestrare e valutare un decoder di posizione
%  del target di reach basato sul firing rate medio, per
%  diverse condizioni sperimentali (free-gaze, motor-only,
%  controlled-gaze, gaze-only).
%
%  Per ciascun file .mat in 'filename':
%    - viene caricato il dataset corrispondente
%    - vengono selezionati automaticamente gli stati PRE e POST
%      (es. 'Pres12', 'Gaze', 'Reach') in base al nome del file
%    - per ogni trial e per tutti i canali/array vengono estratti
%      i conteggi di spike in due finestre temporali:
%         * fase PRE  : ultimi 100 ms dello stato PRE
%         * fase POST : primi 500 ms dello stato POST
%    - le due finestre (PRE+POST) vengono concatenate nel tempo
%      e si calcola il firing rate medio per canale nel trial
%      (Hz), ottenendo la matrice delle feature X
%    - si costruisce il vettore etichette Y (Target_ID)
%    - si esegue una k-fold cross-validation (k = 5) tramite
%      la funzione esterna svm_cv per addestrare un SVM multiclass
%      (fitcecoc con kernel RBF) e stimare l’accuratezza media
%    - si visualizza la confusion matrix per ciascuna condizione
%
%  Al termine:
%    - le accuratezze di tutte le condizioni vengono raccolte
%      in un barplot con indicazione del livello di chance.
%
%  Finestra temporale analizzata:
%     w = [-100, +500] ms implementata come:
%       - ultimi 100 ms dello stato PRE
%       - primi 500 ms dello stato POST
%
%  Prima dell’esecuzione, modificare se necessario:
%    - l’elenco dei file .mat nel cell array 'filename'
%    - i parametri bin_size, period_pre, period_post
%    - la funzione svm_cv (k-fold, tipo di SVM, ecc.)
% =========================================================

clearvars
close all
clc

n_sets = 6;
n_arrays = 2;
n_trials = 32;
bin_size = 0.02;
period_pre = 0.1;
period_post = 0.5;

filename = { ...
    '../00_Data_extraction/free-gaze_BCI02.mat' ...
    '../00_Data_extraction/motor_BCI02.mat' ...
    '../00_Data_extraction/controlled_BCI02.mat' ...
    '../00_Data_extraction/gaze_BCI02.mat' ...
};

%% Costruzione vettore Y e matrice X per SVM + k-fold cross-validation
cm_all  = cell(numel(filename),1);
acc_all = cell(numel(filename),1);
for d = 1:numel(filename) 
    fprintf('\nDataset: %s\n', filename{d}); 
    load(filename{d});

    [~, baseName, ext] = fileparts(filename{d});
    ds_name = [baseName ext];

    if strcmp(ds_name, 'controlled_BCI02.mat')
        PRE  = "Gaze";
        POST = "Reach";
    elseif strcmp(ds_name, 'gaze_BCI02.mat')
        PRE  = "Pres12";
        POST = "Gaze";
    else
        PRE  = "Pres12";
        POST = "Reach";
    end

    Y = [];
    for set = 1:n_sets
        Y = [Y; [data(set).Data(1).Resampled.Target_ID]'];
    end

    idx_pres  = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == PRE);
    idx_reach = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == POST);

    n_bins_pre  = round(period_pre/bin_size);
    n_bins_post = round(period_post/bin_size);

    start_pre = size(data(1).Data(1).Resampled(1).Task_states{idx_pres,2}, 1) - n_bins_pre + 1;
    end_post  = n_bins_post;

    j = 1;
    X = cell(n_trials*n_sets,1);

    for set = 1:n_sets
        for trial = 1:n_trials

            tmp_pre = [];
            tmp_post = [];

            for array = 1:n_arrays
                tmp_pre  = [tmp_pre,  data(set).Data(array).Resampled(trial).Task_states{idx_pres,2}(start_pre:end, :)];
                tmp_post = [tmp_post, data(set).Data(array).Resampled(trial).Task_states{idx_reach,2}(1:end_post, :)];
            end

            matrix = [tmp_pre; tmp_post];
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

    cc = confusionchart(cm, classnames, 'Normalization','row-normalized', ...
        'RowSummary','off', 'ColumnSummary','off');
    cc.Title = sprintf('Confusion Matrix - Location Decoder (Accuracy: %.2f%%)', acc*100);
    cc.XLabel = 'Predicted Target';
    cc.YLabel = 'True Target';

    cm_all{d,1} = cm;
    acc_all{d,1} = acc;

end

%% Barplot
figure('Color','w');

acc = cell2mat(acc_all)*100;  % accuracy in percentuale
nClassi = numel(unique(Y));   % numero di target/classi
chance = (1/nClassi)*100;     % livello di chance

b = bar(acc, 0.5, 'FaceColor', [0.2 0.4 0.7], 'EdgeColor', 'none');
ylim([0 100]);
ylabel('Performance (%)');
xlabel('Condition');
xticks(1:4);
xticklabels({'Free-gaze', 'Motor-only', 'Controlled-gaze', 'Gaze-only'});
title('SVM overall accuracy across conditions');

% Linea di chance level
hold on;
yline(chance, '-k', sprintf('Chance level (%.1f%%)', chance), ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', ...
    'FontSize', 9, 'LineWidth', 1.2);

grid on;
box off;
