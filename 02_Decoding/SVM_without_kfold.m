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

n_sets = 6; 
n_arrays = 2;
n_trials = 32;  
bin_size = 0.02;
period_pre = 0.1; 
period_post = 0.5; 

filename = '../00_Data_extraction/free-gaze_BCI02.mat';
load(filename)

[~, baseName, ext] = fileparts(filename);
ds_name = [baseName ext];

if strcmp(ds_name, 'controlled_BCI02.mat')
        PRE = "Gaze";
        POST = "Reach";
    elseif strcmp(ds_name, 'gaze_BCI02.mat')
        PRE = "Pres12";
        POST = "Gaze";
    else
        PRE = "Pres12";
        POST = "Reach";
end

%% Costruzione vettore Y e matrice X per SVM
Y = []; 
for set = 1:n_sets
    Y = [Y; [data(set).Data(1).Resampled.Target_ID]']; 
end 

idx_pres = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == PRE); 
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
            tmp_pre = [tmp_pre, data(set).Data(array).Resampled(trial).Task_states{idx_pres,2}(start_pre:end, :)]; 
            tmp_post = [tmp_post, data(set).Data(array).Resampled(trial).Task_states{idx_reach,2}(1:end_post, :)]; 
        end 
        matrix = [tmp_pre; tmp_post]; 
        X{j} = mean(matrix./bin_size,1);
        j = j + 1; 
    end   
end 
X = cell2mat(X); 

%% Divisione training/test
cv = cvpartition(Y,'HoldOut',0.3,'Stratify',true);
idxTrain = training(cv); 
idxTest = test(cv);

Xtrain = X(idxTrain, :);   
Xtest = X(idxTest, :);   

Ytrain = Y(idxTrain);
Ytest  = Y(idxTest);

%% Addestramento decoder (SVM multiclass)
t = templateSVM('KernelFunction','rbf','KernelScale','auto','Standardize',true);
Msvm = fitcecoc( ...
    Xtrain, Ytrain, ...
    'Learners', t, ...
    'Coding', 'onevsall');

%% Predizione
Ypred = predict(Msvm, Xtest);

%% Accuratezza
accuracy = mean(Ypred == Ytest);
disp(['Decoder accuracy: ', num2str(accuracy*100, '%.2f'), '%']);

% Confusion matrix con percentuali
classes = unique(Y);
classLabels = arrayfun(@(c) sprintf('Target %d', c), classes, 'UniformOutput', false);
Ytest_cat = categorical(Ytest, classes, classLabels);
Ypred_cat = categorical(Ypred, classes, classLabels);

%% Figure
figure('Color','w');
cm = confusionchart(Ytest_cat, Ypred_cat);
cm.Normalization = 'row-normalized';  % percentuali su ogni riga
cm.Title = sprintf('Confusion Matrix - Location Decoder (Accuracy: %.2f%%)', accuracy*100);
cm.XLabel = 'Predicted Target';
cm.YLabel = 'True Target';

