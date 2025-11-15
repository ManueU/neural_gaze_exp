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

condLab = {'Free-gaze', 'Gaze-on-center', 'Gaze-on-target', 'Gaze-only'};

%% Costruzione vettore Y e matrice X per SVM + k-fold cross-validation
X_all   = cell(numel(filename),1);
Y_all   = cell(numel(filename),1);
classes_all = cell(numel(filename),1);
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
    classes = unique(Y);

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

    % Salvo X e Y per questa condizione
    X_all{d,1} = X;
    Y_all{d,1} = Y;

end

nClassi = numel(classes);
acc_cross = nan(numel(filename), numel(filename) );
cm_cross  = cell(numel(filename), numel(filename) );

t = templateSVM('KernelFunction','rbf', 'KernelScale','auto', 'Standardize',true);
for iTr = 1:numel(filename) 
    fprintf('\n### TRAIN condizione: %s ###\n', condLab{iTr});

    X_train = X_all{iTr};
    Y_train = Y_all{iTr};

    % Addestramento decoder
    Msvm = fitcecoc(X_train, Y_train, ...
        'Learners', t, ...
        'Coding', 'onevsall', ...
        'ClassNames', classes);

    for iTe = 1:numel(filename) 
        fprintf('   -> TEST su: %s\n', condLab{iTe});

        X_test = X_all{iTe};
        Y_test = Y_all{iTe};

        % Predizione
        Y_pred = predict(Msvm, X_test);

        % Confusion matrix ordinata
        cm = confusionmat(Y_test, Y_pred, 'Order', classes);
        cm_cross{iTr, iTe} = cm;

        % Accuracy globale
        acc_cross(iTr, iTe) = sum(diag(cm)) / sum(cm(:));
    end
end

%% --- (3) HEATMAP DEL CROSS-DECODING ---
figure('Color','w');
imagesc(acc_cross * 100);
colormap(copper);  
colorbar;
caxis([0 100]);
axis square;

xticks(1:numel(filename) );
yticks(1:numel(filename) );
xticklabels(condLab);
yticklabels(condLab);

xlabel('Test condition');
ylabel('Train condition');
title('Cross-decoding accuracy (%)');

%% --- (4) Barplot del cross-decoding (colori pastello) ---
figure('Color','w');

% accuracy in percentuale
b = bar(1:numel(filename) , acc_cross * 100, 'grouped', 'EdgeColor', 'none');
ylabel('Accuracy (%)');
xlabel('TRAIN condition');
title('Cross-decoding accuracy');
xticks(1:numel(filename) );
xticklabels(condLab);
ylim([0 100]);

% --- Palette pastello ---
greenPastel = [
    0.05 0.20 0.10   % night green (dark forest)
    0.25 0.45 0.30   % darker soft green
    0.50 0.75 0.55   % medium pastel green
    0.85 0.95 0.85   % very light green
];


% Assegno i colori alle barre (riciclandoli se servono)
for k = 1:numel(filename) 
    colorIndex = mod(k-1, size(greenPastel,1)) + 1;
    b(k).FaceColor = greenPastel(colorIndex, :);
end

% Chance line
chance = (1/nClassi) * 100;
hold on;
yline(chance, '--k', 'Chance', ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','top', ...
    'LineWidth', 0.75);

% Legenda
legend(condLab, ...
    'Location','northoutside', ...
    'Orientation','horizontal', ...
    'Box','off');

grid on;
box off;
clear set
set(gca,'FontSize',11);

%% --- Cross-decoding matrix in stile confusion matrix ---
figure('Color','w');

imagesc(acc_cross * 100);          % percentuale
blueNightMap = [
    1.00 1.00 1.00   % white
    0.85 0.92 1.00   % very light blue
    0.45 0.62 0.80   % steel blue
    0.05 0.10 0.25   % night blue
];
colormap(interp1(1:4, blueNightMap, linspace(1,4,256)));


colorbar;
caxis([0 100]);
axis square;

xticks(1:numel(filename));
yticks(1:numel(filename));
xticklabels(condLab);
yticklabels(condLab);

xlabel('Test condition');
ylabel('Train condition');
title('Cross-decoding accuracy (%)');

set(gca, 'FontSize', 12);

