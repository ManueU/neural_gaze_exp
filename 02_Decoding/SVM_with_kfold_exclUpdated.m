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
n_trials = 16;
bin_size = 0.02;
period_pre = 1.0;
period_post = 0.0;

filename = {
           % '../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat',... 
           % '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat',...
           % '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat'
           % ,...
           '../00_Data_extraction/BCI02_Session_0924/gaze_BCI02_exclUpdated.mat'
           };

%% Costruzione vettore Y e matrice X per SVM + k-fold cross-validation
cm_all  = cell(numel(filename),1);
acc_all = cell(numel(filename),1);
metrics_all = cell(numel(filename),1);

for d = 1:numel(filename) 
    fprintf('\nDataset: %s\n', filename{d}); 
    load(filename{d});

    [~, baseName, ext] = fileparts(filename{d});
    ds_name = [baseName ext];

    if strcmp(ds_name, 'controlled_BCI02_exclUpdated.mat')
        PRE  = "Gaze";
        POST = "Reach";
    elseif strcmp(ds_name, 'gaze_BCI02_exclUpdated.mat')
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

    % Controllo sbilanciamento e k-fold adattivo
    Ycat = categorical(Y);
    counts = countcats(Ycat);
    minCount = min(counts);
    
    k_fold = min(6, minCount);
    n_repeats = 20;
    seed = 0;
    
    if k_fold < 2
        error('Una o più classi hanno <2 campioni: impossibile fare k-fold CV.');
    end
    
    fprintf('Class counts: '); fprintf('%d ', counts); fprintf('\n');
    fprintf('Using repeated %d-fold CV with %d repeats\n', k_fold, n_repeats);
    
    % [acc, cm, metrics] = svm_cv(X, Y, k_fold);
    [results, cm_pooled] = svm_repeated_cv(X, Y, k_fold, n_repeats, seed);
    fprintf('Accuracy: %.3f ± %.3f\n', mean(results.acc), std(results.acc));
    fprintf('Balanced accuracy: %.3f ± %.3f\n', mean(results.balAcc), std(results.balAcc));
    fprintf('Macro-F1: %.3f ± %.3f\n', mean(results.macroF1), std(results.macroF1));
    
    classes = categories(categorical(Y));

    figure('Color','w');
    cc = confusionchart(cm_pooled, classes, ...
        'Normalization','row-normalized', ...
        'RowSummary','off', ...
        'ColumnSummary','off');

    cc.GridVisible = 'off';
    cc.DiagonalColor = [0 0.35 0.7];
    cc.OffDiagonalColor = [0.8 0.8 0.8];
    cc.XLabel = '\bf Predicted Target';
    cc.YLabel = '\bf True Target';

    cm_all{d} = cm_pooled;
    acc_all{d} = results.acc;
    metrics_all{d} = results;

end

%% Figure (1) - Balanced accuracy per condizione + recall per classe
% ============================================================
% Questo grafico riassume la performance del decoder SVM nelle
% diverse condizioni sperimentali usando come metrica principale
% la balanced accuracy.
%
% Per ciascuna condizione:
%   - la barra rappresenta la balanced accuracy media ottenuta
%     sui repeat della cross-validation;
%   - la error bar rappresenta la deviazione standard tra i repeat,
%     cioè la variabilità della performance dovuta alle diverse
%     partizioni train/test;
%   - i puntini mostrano il recall delle singole classi (target),
%     calcolato dalla confusion matrix pooled, e permettono di
%     visualizzare se il decoding è omogeneo tra i diversi target
%     oppure se alcune classi sono sistematicamente più facili o
%     più difficili da classificare.
%
% La linea tratteggiata indica il chance level teorico, pari a
% 1 / numero di classi.
%
% Questa figura permette quindi di confrontare sia la performance
% media globale tra condizioni sia la distribuzione della
% performance sulle singole classi.
% ============================================================

nCond   = numel(cm_all);
nClassi = size(cm_all{1}, 1);

balAcc_mean = zeros(nCond,1);
balAcc_std  = zeros(nCond,1);

classAcc = zeros(nCond,nClassi);

for d = 1:nCond
    
    % balanced accuracy dai repeat
    balAcc = metrics_all{d}.balAcc;
    
    balAcc_mean(d) = mean(balAcc) * 100;
    balAcc_std(d)  = std(balAcc) * 100;
    
    % recall per classe dalla confusion matrix pooled
    cm = cm_all{d};
    classAcc(d,:) = diag(cm) ./ sum(cm,2);
    
end

figure('Color','w')

bar(1:nCond, balAcc_mean, 0.5, ...
    'FaceColor',[0.68 0.41 0.44], ...
    'EdgeColor','none')

hold on

errorbar(1:nCond, balAcc_mean, balAcc_std, ...
    'k','LineStyle','none','LineWidth',1.5)

% chance level
chance = (1/nClassi) * 100;
yline(chance,'--k','Chance', ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','top', ...
    'LineWidth',0.5)

% punti delle classi
for d = 1:nCond
    
    x_vals = d + (rand(1,nClassi)-0.5)*0.15;
    y_vals = classAcc(d,:) * 100;
    
    scatter(x_vals,y_vals,30,'filled',...
        'MarkerFaceColor',[0.90 0.65 0.70],...
        'MarkerEdgeColor','none',...
        'MarkerFaceAlpha',0.6)
end

ylim([0 100])
ylabel('Balanced accuracy (%)')
xlabel('Condition')

xticks(1:nCond)
xticklabels({'Free-gaze','Gaze-on-center','Gaze-on-target','Gaze-only'})

title('SVM balanced accuracy across conditions')

grid on
box off


%% Figure (2) - Violin plot balanced accuracy across conditions
% ============================================================
% Questo grafico mostra la distribuzione completa della balanced
% accuracy nelle diverse condizioni sperimentali, usando i valori
% ottenuti dai repeat della cross-validation.
%
% Per ciascuna condizione, il violin plot rappresenta la densità
% dei valori di balanced accuracy:
%   - violini più larghi indicano una maggiore concentrazione di
%     repeat in quel range di performance;
%   - la forma del violino permette di visualizzare variabilità,
%     asimmetrie o eventuali valori estremi.
%
% A differenza del barplot, questa figura non riassume soltanto
% la media, ma mostra l’intera distribuzione della performance
% del decoder, offrendo una rappresentazione più informativa
% della stabilità del decoding.
%
% La linea tratteggiata indica il chance level teorico, pari a
% 1 / numero di classi.
%
% Questa figura è utile per confrontare non solo il livello medio
% di performance tra condizioni, ma anche la loro variabilità e
% robustezza rispetto alla suddivisione dei dati nei fold.
% ============================================================

nCond = numel(metrics_all);

balAcc_all = [];
group = [];

for d = 1:nCond
    
    balAcc = metrics_all{d}.balAcc * 100;
    
    balAcc_all = [balAcc_all; balAcc];
    group = [group; d*ones(numel(balAcc),1)];
    
end

figure('Color','w')

violinplot(balAcc_all, group)

ylim([0 100])
ylabel('Balanced accuracy (%)')
xlabel('Condition')

xticklabels({'Free-gaze','Gaze-on-center','Gaze-on-target','Gaze-only'})

hold on

% chance level
nClassi = size(cm_all{1},1);
chance = (1/nClassi) * 100;

yline(chance,'--k','Chance', ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','top', ...
    'LineWidth',0.5)

title('Distribution of decoding performance (repeated CV)')

grid on
box off

%% Figure (3) - Representational Dissimilarity Matrix (RDM)
% ============================================================
% Questa matrice quantifica quanto i pattern di errore del decoder
% siano simili o diversi tra le condizioni sperimentali.
%
% L’idea è confrontare non tanto la performance media del decoder,
% ma la struttura delle sue confusioni: due condizioni risultano
% simili se il decoder tende a confondere le stesse classi in modo
% analogo.
%
% Procedura:
%   1) ogni confusion matrix viene normalizzata per riga, così che
%      ogni classe contribuisca in modo comparabile;
%   2) la diagonale principale viene rimossa, per escludere le
%      predizioni corrette e considerare solo gli errori;
%   3) ciascuna matrice viene convertita in un vettore;
%   4) si calcola la correlazione di Spearman tra i vettori delle
%      diverse condizioni;
%   5) la dissimilarità è definita come:
%
%           RDM = 1 - correlazione
%
% Interpretazione dei valori:
%   - valori vicini a 0 indicano pattern di errore molto simili;
%   - valori più alti indicano pattern di errore più differenti.
%
% Questa figura permette quindi di valutare se condizioni diverse
% condividano una struttura simile di errori del decoder, suggerendo
% una possibile somiglianza nella rappresentazione neurale dei target.
% ============================================================

nCond = numel(cm_all);
vec_all = cell(nCond,1);

for d = 1:nCond
    
    cm = cm_all{d};
    
    % normalizzazione per riga
    row_sum = sum(cm,2);
    cm_norm = cm ./ max(row_sum,1);
    
    % rimuove la diagonale (corretti)
    cm_norm(1:size(cm_norm,1)+1:end) = 0;
    
    % vettorizzazione
    vec_all{d} = cm_norm(:);
    
end

% correlazioni
R = zeros(nCond);

for i = 1:nCond
    for j = 1:nCond
        
        R(i,j) = corr(vec_all{i}, vec_all{j}, ...
                      'type','Spearman','rows','complete');
        
    end
end

% dissimilarità
RDM = 1 - R;

% Plot

figure('Color','w')

imagesc(RDM)

colormap(flipud(gray))
colorbar

caxis([0 1])

axis square

xticks(1:nCond)
yticks(1:nCond)

cond_names = {'Free-gaze','Gaze-on-center','Gaze-on-target','Gaze-only'};

xticklabels(cond_names)
yticklabels(cond_names)

title('Representational dissimilarity matrix (decoder error patterns)')