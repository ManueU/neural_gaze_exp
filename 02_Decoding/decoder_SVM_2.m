% =========================================================
% DESCRIZIONE:
%  Script per addestrare e valutare un decoder di posizione
%  del target di reach basato sul firing rate medio.
%
%  Passaggi:
%    - Per ogni trial e per tutti i channels viene estratto il 
%      firing rate in una finestra temporale che unisce:
%         * fase PRE  : ultimi 100 ms prima dell'onset del reach
%         * fase POST : primi 500 ms dopo l'onset del reach
%    - si costruisce così un vettore di feature medio per trial
%    - si costruisce il vettore etichette Y (Target_ID)
%    - si divide il dataset in training/test con cvpartition
%    - si addestra un SVM multiclass (fitcecoc, kernel RBF)
%    - si calcola l'accuratezza sul test set
%    - si visualizza la confusion matrix (normalizzata per riga)
%
%  La finestra temporale analizzata è fissata a:
%     w = [-100, +500] ms rispetto all’onset del reach.
%
%  Prima dell’esecuzione, modificare se necessario:
%    - filename : nome del file .mat da caricare
%    - PRE, POST: etichette degli stati di interesse
%                 (ad es. 'Pres12', 'Gaze', 'Reach')
% =========================================================

clearvars
close all
clc
load('controlled_BCI02.mat')

n_sets   = 6; 
n_trials = 32*ones(1,n_sets);  
bin_size = 0.02;

%Stati
NAME_PRE     = "Gaze"; 
NAME_REACH   = "Reach"; 
period_pre   = 0.1; 
period_reach = 0.5; 

% Preparatory window
idx_pres = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == NAME_PRE); 
start_pres = size(data(1).Data(2).Resampled(1).Task_states{idx_pres,2}, 1) - (period_pre/bin_size); 

% Movement window
idx_reach = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == NAME_REACH); 
end_reach = period_reach/bin_size;

% Labels
tmp = []; 
for set = 1:n_sets
    Y = [tmp; [data(set).Data(1).Resampled.Target_ID]']; 
    tmp = Y;
end 

% Data
j = 1; 
X = cell(sum(n_trials),1);
for set = 1:n_sets
    for trial = 1:n_trials(set)
        tmp_pres = []; 
        tmp_reach = []; 
        for array = 1:2
            tmp_pres = [tmp_pres, data(set).Data(array).Resampled(trial).Task_states{idx_pres,2}(start_pres:end, :)]; 
            tmp_reach = [tmp_reach, data(set).Data(array).Resampled(trial).Task_states{idx_reach,2}(1:end_reach, :)]; 
        end 
        matrix = [tmp_pres; tmp_reach]; 
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
cc = confusionchart(cm, classnames, 'Normalization','row-normalized', 'RowSummary','off', 'ColumnSummary','off');
cc.Title  = sprintf('Confusion Matrix - Location Decoder (Accuracy: %.2f%%)', acc*100);
cc.XLabel = 'Predicted Target';
cc.YLabel = 'True Target';

