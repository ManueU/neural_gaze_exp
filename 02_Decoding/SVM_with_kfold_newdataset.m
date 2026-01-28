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
%    - si esegue una k-fold cross-validation (k = n) tramite
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
% close all
clc

n_sets = 6;
n_arrays = 2;
n_trials = 32;
bin_size = 0.02;
period_pre = 0.1;
period_post = 0.5;

% filename = { ...
%     '../00_Data_extraction/free-gaze_BCI02.mat' ...
%     '../00_Data_extraction/motor_BCI02.mat' ...
%     '../00_Data_extraction/controlled_BCI02.mat' ...
%     '../00_Data_extraction/gaze_BCI02.mat' ...
% };

filename = {'../00_Data_extraction/free-gaze_BCI02_withtracker_exclUpdated.mat',... 
                   '../00_Data_extraction/motor_BCI02_withtracker_exclUpdated.mat',...
                   '../00_Data_extraction/controlled_BCI02_withtracker_exclUpdated.mat',...
                   '../00_Data_extraction/gaze_BCI02_withtracker_exclUpdated.mat'};

%% Costruzione vettore Y e matrice X per SVM + k-fold cross-validation
cm_all  = cell(numel(filename),1);
acc_all = cell(numel(filename),1);
for d = 1:numel(filename) 
    fprintf('\nDataset: %s\n', filename{d}); 
    load(filename{d});

    [~, baseName, ext] = fileparts(filename{d});
    ds_name = [baseName ext];

    if strcmp(ds_name, 'controlled_BCI02_withtracker_exclUpdated.mat')
        PRE  = "Gaze";
        POST = "Reach";
    elseif strcmp(ds_name, 'gaze_BCI02_withtracker_exclUpdated.mat')
        PRE  = "Pres12";
        POST = "Gaze";
    else
        PRE  = "Pres12";
        POST = "Reach";
    end

    Y = [];
    for set = 1:n_sets
        idx = [data(set).Data(1).Interp.Excluded] == 0;
        Y = [Y; [data(set).Data(1).Interp(idx).Target_ID]'];
    end

    idx_pres  = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == PRE);
    idx_reach = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == POST);

    n_bins_pre  = round(period_pre/bin_size);
    n_bins_post = round(period_post/bin_size);

    start_pre = size(data(1).Data(1).Interp(1).Task_states{idx_pres,2}, 1) - n_bins_pre + 1;
    end_post  = n_bins_post;
    n_valid = sum( arrayfun(@(s) sum([data(s).Data(1).Interp.Excluded] == 0), 1:n_sets) );


    j = 1;
    X = cell(n_valid, 1);

    for set = 1:n_sets
        for trial = 1:n_trials

            tmp_pre = [];
            tmp_post = [];
            if data(set).Data(1).Interp(trial).Excluded == 0
                for array = 1:n_arrays
                    tmp_pre  = [tmp_pre,  data(set).Data(array).Interp(trial).Task_states{idx_pres,2}(start_pre:end, :)];
                    tmp_post = [tmp_post, data(set).Data(array).Interp(trial).Task_states{idx_reach,2}(1:end_post, :)];
                end

                matrix = [tmp_pre; tmp_post];
                X{j} = mean(matrix./bin_size,1);
                j = j + 1;
            end
        end
    end

    X = cell2mat(X);

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
    [acc, cm, metrics] = svm_cv(X, Y, k_fold);
    
    figure('Color','w');
    classes = unique(Y);
    classnames = arrayfun(@(c) sprintf('Target %d', c), classes, 'UniformOutput', false);

    cc = confusionchart(cm, classes, 'Normalization','row-normalized', ...
        'RowSummary','off', 'ColumnSummary','off');
    cc.GridVisible = 'off';
    cc.DiagonalColor    = [0 0.35 0.7];
    cc.OffDiagonalColor = [0.8 0.8 0.8];    % grigio chiaro
    cc.Title = sprintf('Acc: %.2f%% | BalAcc: %.2f%% | MacroF1: %.2f%%', ...
                   acc*100, metrics.balancedAccuracy*100, metrics.macroF1*100);
    cc.XLabel = '\bf Predicted Target';
    cc.YLabel = '\bf True Target';

    cm_all{d,1} = cm;
    acc_all{d,1} = acc;

end

%% Figure (1) - Barplot overall accuracy per ciascuna condizione
% ============================================================
% Questo grafico mostra la performance globale del decoder SVM
% nelle 4 condizioni sperimentali. Per ciascun dataset viene
% riportata l’accuracy totale sul test set, cioè:
%
%       accuracy = (# predizioni corrette) / (totale)
% 
% È una misura generale, ma non dettagliata.
% ============================================================

figure('Color','w');

acc_overall = cell2mat(acc_all) * 100;   % accuracy in percentuale
nCond   = numel(acc_overall);            % numero condizioni (4)
nClassi = size(cm_all{1}, 1);            % numero di classi (dalle confusion matrix)
chance  = (1/nClassi) * 100;             % livello di chance in %

b = bar(1:nCond, acc_overall, 0.5, ...
        'FaceColor', [0.68, 0.41, 0.44], ...
        'EdgeColor', 'none'); 
ylim([0 100]);
ylabel('Performance (%)');
xlabel('Condition');
xticks(1:nCond);
xticklabels({'Free-gaze', 'Gaze-on-center', 'Gaze-on-target', 'Gaze-only'});
title('SVM overall accuracy across conditions');

hold on;
yline(chance, '--k', 'Chance','LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LineWidth', 0.5);

grid on;
box off;


%% Figure (2) - Balanced accuracy delle singole class
% ============================================================
% La balanced accuracy è la media dell'accuratezza per ciascuna
% classe (target), quindi:
%
%       balAcc = mean( diag(cm) ./ sum(cm,2) )
%
% Questo evita che classi sbilanciate dominino la
% performance totale.
%
% Sul barplot vengono sovrapposti i punti relativi all’accuracy
% di ogni singola classe, così da visualizzare come varia la
% performance tra i diversi target all'interno di ogni condizione.
%
% La figura ci dice se la performance è omogenea su tutti i target o no.
% ============================================================

nCond   = numel(cm_all);
nClassi = size(cm_all{1}, 1);
classAcc = zeros(nCond, nClassi);  
for d = 1:nCond
    cm = cm_all{d};                          % confusion matrix della condizione d
    perClass = diag(cm) ./ sum(cm, 2);       % accuracy per classe (riga normalizzata)
    classAcc(d, :) = perClass.';             % 1 x nClassi
end
balancedAcc = mean(classAcc, 2) * 100;       

figure('Color','w');
b = bar(1:nCond, balancedAcc, 0.5, ...
        'FaceColor', [0.68, 0.41, 0.44], ...
        'EdgeColor', 'none');
hold on;

% Chance level
chance = (1/nClassi) * 100;
yline(chance, '--k', 'Chance','LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'LineWidth', 0.5);

% Punti delle singole classi sopra ogni barra
for d = 1:nCond
    x_vals = d * ones(1, nClassi); 
    x_vals = x_vals + (rand(1, nClassi)-0.5)*0.15;

    y_vals = classAcc(d, :) * 100;        
    
    scatter(x_vals, y_vals, 30, 'filled', ...
            'MarkerFaceColor', [0.90, 0.65, 0.70], ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6);
end

ylim([0 100]);
ylabel('Performance (%)');
xlabel('Condition');
xticks(1:nCond);
xticklabels({'Free-gaze', 'Gaze-on-center', 'Gaze-on-target', 'Gaze-only'});
title('SVM balanced accuracy across conditions');
grid on;
box off;

%% Figure (3) - Representational Dissimilarity Matrix (RDM)
% ============================================================
% Questa matrice mostra la dissimilarità tra le condizioni
% sperimentali sulla base dei pattern di confusione del decoder.
% Risponde alla domanda: Il decoder sbaglia nello stesso modo 
% in condizioni diverse?
%
% Procedura:
%   1) Ogni confusion matrix viene normalizzata per riga.
%   2) Appiattita in un vettore (pattern degli errori).
%   3) Si calcola la correlazione tra vettori di condizioni diverse.
%   4) La dissimilarità viene ottenuta come: RDM = 1 - correlazione.
%
% Valori:
%   - 0   = pattern molto simili (neri)
%   - 1   = pattern diversi
%   - >1  = anticorrelazione (raramente)
%
% Questo grafico evidenzia quali condizioni condividono una simile
% struttura degli errori del decoder, riflettendo similarità nelle
% rappresentazioni neurali dei target.
% ============================================================

nCond = numel(cm_all);
vec_all = cell(nCond,1);
for d = 1:nCond
    cm = cm_all{d};
    cm_norm = cm ./ sum(cm,2);  
    vec_all{d} = cm_norm(:);   
end

R = zeros(nCond);
for i = 1:nCond
    for j = 1:nCond
        R(i,j) = corr(vec_all{i}, vec_all{j}, 'type','Pearson');
    end
end
RDM = 1 - R;

figure('Color','White');
imagesc(RDM);
colormap(gray); colorbar;
title('Representational dissimilarity (RDM)');
xticks(1:4); yticks(1:4);
xticklabels({'Free-gaze', 'Gaze-on-center', 'Gaze-on-target', 'Gaze-only'});
yticklabels({'Free-gaze', 'Gaze-on-center', 'Gaze-on-target', 'Gaze-only'});
