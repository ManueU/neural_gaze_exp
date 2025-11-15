% =========================================================
% DESCRIZIONE:
%  Script per eseguire una decodifica temporale del Target_ID
%  tramite SVM multiclass, utilizzando il firing rate medio
%  estratto da finestre temporali scorrevoli su dati neurali.
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
%    - si costruisce il vettore delle etichette Y (Target_ID),
%      concatenando i Target_ID di tutti i set
%    - si suddivide il dataset in training (70%) e test (30%)
%      tramite cvpartition con HoldOut stratificato
%    - si addestra, per la finestra corrente, un classificatore
%      SVM multiclass (fitcecoc + templateSVM RBF)
%    - si stimano:
%         * l’accuratezza complessiva sul set di test
%         * l’accuratezza per ciascuna classe (diag normalizzata
%           della confusion matrix)
%
%  Visualizzazione finale:
%    - andamento temporale dell’accuratezza complessiva e delle
%      accuratezze per classe (curve smussate nel tempo con
%      smoothdata gaussiano)
%    - indicazione degli eventi/stati del task tramite linee
%      verticali e relative etichette
%    - indicazione del livello di chance (1 / n_classes)
%    - asse temporale basato sul centro di ciascuna finestra
%      scorrevole
%
%  Parametri principali modificabili:
%    - n_sets, n_arrays, n_trials
%    - bin_size, w_length, overlap
%    - proporzione train/test in cvpartition
%    - parametri dell’SVM (KernelFunction, KernelScale, ecc.)
% =========================================================


clearvars
close all
clc

n_sets = 6;
n_arrays = 2;
n_trials = 32;
bin_size = 0.02;

filename = '../00_Data_extraction/controlled_BCI02.mat';
load(filename)

%% Decoding over time with SVM
% bins per trial
N = size(data(1).Data(1).Resampled(1).Trial, 1); 
rec_duration = N*bin_size; 

% finestra scorrevole
w_length = 0.6; 
overlap = 0.5*w_length;
N_w = round(w_length/bin_size);
N_o = round(overlap/bin_size);
n_acc = floor((N - N_w)/(N_w-N_o)) + 1;

% etichette Y
Y = []; 
for set = 1:n_sets
    Y = [Y; [data(set).Data(1).Resampled.Target_ID]']; 
end 
classes = unique(Y,'stable');
n_classes = numel(classes);
classLabels = arrayfun(@(c) sprintf('Target %d', c), classes, 'UniformOutput', false);

rng(0);                       
cv = cvpartition(Y,'HoldOut',0.3,'Stratify',true);
idxTrain = training(cv); 
idxTest = test(cv);
t = templateSVM('KernelFunction','rbf','KernelScale','auto','Standardize',true);    

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
    
    % Divisione training/test
    Xtrain = X(idxTrain, :);   
    Xtest = X(idxTest, :);   
    
    Ytrain = Y(idxTrain);
    Ytest  = Y(idxTest);
    
    % Addestramento decoder (SVM multiclass)
    Msvm = fitcecoc( ...
        Xtrain, Ytrain, ...
        'Learners', t, ...
        'Coding', 'onevsall');
    
    % Predizione
    Ypred = predict(Msvm, Xtest);
    
    % Overall accuracy
    acc_overall{w} = mean(Ypred == Ytest)*100;

    % Accuratezza
    [cm, order] = confusionmat(Ytest, Ypred);
    acc{w} = diag(cm) ./ sum(cm,2) * 100;

    start_w = start_w + (N_w - N_o); 
    end_w = start_w + N_w - 1; 
end 

%% Figure
figure('Color', 'White')
events_time_tmp = []; 
for i = 1:length(data(1).Data(2).Resampled(1).Task_states)
    events_time = [events_time_tmp; size(data(1).Data(2).Resampled(1).Task_states{i,2},1)*bin_size];
    events_time_tmp = events_time; 
end 
increment_times = cumsum(events_time); 
labels = string(data(1).Data(2).Resampled(1).Task_states(:,1));

t = (((0:n_acc-1)*(N_w - N_o)) + N_w/2)*bin_size;
n_targets = numel(order);
acc_mat = zeros(n_acc, n_targets);
for w = 1:n_acc
    acc_mat(w, :) = acc{w}(:).';
end

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
w_smooth = 5;
for c = 1:n_targets
    acc_smooth = smoothdata(acc_mat(:, c), 'gaussian', w_smooth);
    plot(t, acc_smooth, 'LineWidth', 1.0, 'Color', [colors(c,:),  alpha]), hold on
end
acc_smooth_overall = smoothdata(cell2mat(acc_overall), 'gaussian', w_smooth);
plot(t, acc_smooth_overall, 'LineWidth', 1.5, 'Color', 'k'), hold on

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
yline((1/n_classes)*100,'-', 'Chance'); 
xlabel('Time (s)');
ylabel('Accuracy (%)');
xlim([0 rec_duration]);
