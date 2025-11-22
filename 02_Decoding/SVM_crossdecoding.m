% =========================================================
% DESCRIZIONE:
% Script per estrarre feature di firing rate medio e
% valutare il cross-decoding della posizione del target
% di reach tra diverse condizioni sperimentali
% (free-gaze, motor-only, controlled-gaze, gaze-only).
%
% Per ciascun file .mat in 'filename':
% - viene caricato il dataset corrispondente
% - vengono selezionati automaticamente gli stati PRE e POST
% (es. 'Pres12', 'Gaze', 'Reach') in base al nome del file
% - per ogni trial, set e per tutti i canali/array vengono estratti
% i conteggi di spike in due finestre temporali:
% * fase PRE : ultimi 'period_pre' secondi dello stato PRE
% * fase POST : primi 'period_post' secondi dello stato POST
% - le due finestre (PRE+POST) vengono concatenate nel tempo
% e si calcola il firing rate medio per canale nel trial (Hz),
% ottenendo la matrice delle feature X
% - si costruisce il vettore etichette Y (Target_ID) e il vettore
% delle classi distinte 'classes'
%
% Decoding:
% - per ogni condizione di training (riga di 'condLab'):
% * si addestra un SVM multiclass (fitcecoc con kernel RBF)
% sui dati di quella condizione
% * il decoder viene testato su TUTTE le condizioni
% (inclusa quella di training), ottenendo per ogni coppia
% train–test:
% · la confusion matrix (cm_cross)
% · l’accuratezza globale (acc_cross)
%
% Visualizzazioni finali:
% (1) Barplot dell’accuratezza di cross-decoding
% (righe: condizione di training, barre: condizione di test)
% con indicazione del livello di chance (1 / nClassi).
% (2) Barplot delle accuratezze bidirezionali per tutte
% le coppie di condizioni (train A→test B vs train B→test A).
% (3) Matrice di cross-decoding con train sulle righe e test 
% sulle colonne (confusion matrix). 
%
% Finestra temporale analizzata:
% w = [-100, +500] ms 
%
% Prima dell’esecuzione, modificare se necessario:
% - l’elenco dei file .mat nel cell array 'filename'
% - i parametri n_sets, n_arrays, n_trials
% - i parametri bin_size, period_pre, period_post
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

filename = { ...
    '../00_Data_extraction/free-gaze_BCI02.mat' ...
    '../00_Data_extraction/motor_BCI02.mat' ...
    '../00_Data_extraction/controlled_BCI02.mat' ...
    '../00_Data_extraction/gaze_BCI02.mat' ...
};

condLab = {'Free-gaze', 'Gaze-on-center', 'Gaze-on-target', 'Gaze-only'};
nCond = numel(filename);

%% Costruzione vettore Y e matrice X per SVM
X_all   = cell(nCond,1);
Y_all   = cell(nCond,1);
for d = 1:nCond
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
        Y = [Y; [data(set).Data(1).Interp.Target_ID]'];
    end
    classes = unique(Y);

    idx_pres  = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == PRE);
    idx_reach = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == POST);

    n_bins_pre  = round(period_pre/bin_size);
    n_bins_post = round(period_post/bin_size);

    start_pre = size(data(1).Data(1).Interp(1).Task_states{idx_pres,2}, 1) - n_bins_pre + 1;
    end_post  = n_bins_post;

    j = 1;
    X = cell(n_trials*n_sets,1);
    for set = 1:n_sets
        for trial = 1:n_trials

            tmp_pre = [];
            tmp_post = [];

            for array = 1:n_arrays
                tmp_pre  = [tmp_pre,  data(set).Data(array).Interp(trial).Task_states{idx_pres,2}(start_pre:end, :)];
                tmp_post = [tmp_post, data(set).Data(array).Interp(trial).Task_states{idx_reach,2}(1:end_post, :)];
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
acc_cross = nan(nCond, nCond);
cm_cross  = cell(nCond, nCond);

t = templateSVM('KernelFunction','rbf', 'KernelScale','auto', 'Standardize',true);
for iTr = 1:nCond
    fprintf('\n### TRAIN condizione: %s ###\n', condLab{iTr});

    X_train = X_all{iTr};
    Y_train = Y_all{iTr};

    % Addestramento decoder
    Msvm = fitcecoc(X_train, Y_train, ...
        'Learners', t, ...
        'Coding', 'onevsall', ...
        'ClassNames', classes);

    for iTe = 1:nCond
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

%% Figure (1) - Barplot del cross-decoding
figure('Color','w');

% accuracy in percentuale
b = bar(1:nCond, acc_cross * 100, 'grouped', 'EdgeColor', 'none');
ylabel('Accuracy (%)');
xlabel('Train condition');
title('Cross-decoding analysis');
xticks(1:nCond);
xticklabels(condLab);
ylim([0 100]);

greenPastel = [
    0.05 0.20 0.10   % night green (dark forest)
    0.25 0.45 0.30   % darker soft green
    0.50 0.75 0.55   % medium pastel green
    0.85 0.95 0.85   % very light green
];


for k = 1:nCond
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
    'Location','southoutside', ...
    'Orientation','horizontal', ...
    'Box','off');

grid on;
box off;

%% Figure (2)
pairs   = nchoosek(1:nCond, 2);    % ogni riga: [i j]
nPairs  = size(pairs, 1);

vals    = zeros(nPairs, 2);        % colonna 1: i→j, colonna 2: j→i
xLabels = cell(nPairs, 1);

for p = 1:nPairs
    i = pairs(p,1);
    j = pairs(p,2);

    % Accuracy (in %) nelle due direzioni
    vals(p,1) = acc_cross(i,j) * 100;   % train i → test j
    vals(p,2) = acc_cross(j,i) * 100;   % train j → test i

end

figure('Color','w');

b = bar(1:nPairs, vals, 'grouped', 'EdgeColor','none');
ylabel('Accuracy (%)');
title('Cross-decoding analysis');

ax = gca;
ax.XTick = 1:nPairs;
ax.XTickLabel = {};
ylim([0 100]);
yl = ylim;
yText = yl(1)-2;

xx = 1:nPairs;
for p = 1:nPairs
    i = pairs(p,1);
    j = pairs(p,2);

    labelStr = sprintf('%s\nvs.\n%s', condLab{i}, condLab{j});

    text(xx(p), yText, labelStr, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize', 9);
end


bordeaux = [0.45 0.00 0.10];
grigio   = [0.70 0.70 0.70];
b(1).FaceColor = bordeaux;   % train i → test j
b(2).FaceColor = grigio;     % train j → test i

% Linea di chance
if exist('nClassi','var')
    chance = (1/nClassi) * 100;
    hold on;
    yline(chance, '--k', 'Chance', ...
        'LabelHorizontalAlignment','right', ...
        'LabelVerticalAlignment','top', ...
        'LineWidth', 0.75);
end

legend({'Train i → Test j', 'Train j → Test i'}, ...
       'Location','northoutside', ...
       'Orientation','horizontal', ...
       'Box','off');

grid on;
box off;

%% Figure (3) - Cross-decoding matrix (confusion matrix)
figure('Color','w');

imagesc(acc_cross * 100);          
blueNightMap = [
    1.00 1.00 1.00   % white
    0.85 0.92 1.00   % very light blue
    0.45 0.62 0.80   % steel blue
    0.05 0.10 0.25   % night blue
];
colormap(interp1(1:4, blueNightMap, linspace(1,4,256)));


colorbar;
clim([0 100]);
axis square;

xticks(1:nCond);
yticks(1:nCond);
xticklabels(condLab);
yticklabels(condLab);

xlabel('Test condition');
ylabel('Train condition');
title('Cross-decoding accuracy (%)');


