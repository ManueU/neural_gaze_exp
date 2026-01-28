% =========================================================
% DESCRIZIONE:
% Lo script esegue un’analisi dPCA (demixed PCA) su dati neurali
% provenienti da più condizioni sperimentali ("free-gaze", "motor",
% "controlled"), allo scopo di:
%
% - valutare la separabilità dei target di movimento e dello stato di gaze 
%    (condizione comportamentale);
%  - identificare componenti indipendenti per target, condizione, tempo e 
%    interazioni target×condizione;
%  - visualizzare traiettorie neurali proiettate nei sottospazi 
%    Goal-like, Condition-like, Time-dependent e Gaze-dependent.
%
% Il workflow principale comprende:
%   1) Caricamento dei dataset e definizione delle finestre temporali
%      PRE e POST per ciascuna condizione (es. Pres12, Gaze, Reach).
%   2) Estrazione dei firing rate bin-aggregati:
%      - segmentazione delle prove nelle finestre temporali PRE/POST;
%      - media sui trial per ogni target, canale e condizione;
%      - costruzione della matrice completa [tempo × target × neuroni].
%   3) Filtraggio dei neuroni poco informativi tramite:
%      - frequenza di firing media minima;
%      - varianza minima della risposta neurale.
%   4) Costruzione del tensore X per dPCA con struttura:
%         [neuroni × tempo × target × condizione]
%      e conversione nel formato richiesto da dPCA:
%         firingRatesAverage = [neuroni × target × condizione × tempo]
%   5) Esecuzione della dPCA senza regularizzazione:
%      - calcolo di W, V e whichMarg;
%      - quantificazione della varianza spiegata per ciascuna marginalizzazione:
%           Target, Condition (Gaze), Time, Interazione T×C.
%   6) Estrazione delle componenti latenti Z mediante proiezione W·X:
%      - ricostruzione degli andamenti temporali delle componenti;
%      - selezione automatica delle componenti più informative per
%        sottospazi Goal e Gaze (TDR).
%   7) Generazione di figure:
%       • Traiettorie 3D nel sottospazio Goal/Gaze.
%       • Evoluzione temporale delle componenti Time-, Target- e
%         Gaze-dependent.
%       • Istogramma della varianza spiegata da Target, Gaze e Tempo.
%       • Analisi del contributo di ciascuna condizione per target (e viceversa).
%
% Prima dell’esecuzione, modificare i parametri:
%   filename      : lista dei dataset da caricare (.mat)
%   sets          : set/partizioni delle prove da includere
%   period_pre    : durata finestra PRE (s)
%   period_post   : durata finestra POST (s)
%   bin_size      : ampiezza dei bin (s)
%   minMeanFR     : firing rate minimo per mantenere un neurone
%   minVarFR      : varianza minima per mantenere un neurone
% =========================================================

% 1: medial arm 
% 2: lateral hand 

clearvars -except mean_baseline std_baseline
% close all
clc

sets = [2,4,5,6]; 
n_sets = numel(sets); 
n_arrays = 2;          
n_channels = 96;    
n_targets = 8; 
n_trials = 32; 
bin_size = 0.02;       
period_pre = 1.0;     % intervallo migliore 
period_post = 0.5;        

filename = { ...
    '../00_Data_extraction/free-gaze_BCI02.mat', ...
    '../00_Data_extraction/motor_BCI02.mat',...
    '../00_Data_extraction/controlled_BCI02.mat'};  


%% Costruzione matrice condizione (media su trial, finestra PRE/POST)

condition = []; 
condition_sep = cell(numel(filename),1);
for d = 1:numel(filename) 
    fprintf('\nDataset: %s\n', filename{d}); 
    load(filename{d});   

    [~, baseName, ext] = fileparts(filename{d});
    ds_name = [baseName ext];

    % nomi stati PRE / POST a seconda del file
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

    idx_pres  = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == PRE); 
    idx_reach = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == POST); 
    if isempty(idx_pres) || isempty(idx_reach)
        error('Stati PRE/POST non trovati: controlla PRE/POST');
    end
    
    % numero di bin PRE/POST
    sigma_bins = 25 / 20;   % = 1.25

    pre_bins  = max(1, round(period_pre/bin_size));
    post_bins = max(1, round(period_post/bin_size));
    
    pca_matrix_array = []; 
    for array = 1:n_arrays
        pca_matrix = []; 
    
        for channel = 1:n_channels
            zscored_by_targets = []; 
            ch_global = (array-1)*n_channels + channel;

            for target = 1:n_targets
                firing_rate = [];
    
                for set_ = 1:n_sets
                    set = sets(set_);
                    idx = find([data(set).Data(array).Interp.Target_ID] == target); 
    
                    for j = 1:length(idx)
                        tmp_pre  = data(set).Data(array).Interp(idx(j)).Task_states{idx_pres, 2}(end-pre_bins+1:end, channel); 
                        tmp_post = data(set).Data(array).Interp(idx(j)).Task_states{idx_reach,2}(1:post_bins, channel); 
                        % tmp_pre  = data(set).Data(array).Interp(idx(j)).Task_states{idx_pres, 2}(1:end, channel); 
                        % tmp_post = data(set).Data(array).Interp(idx(j)).Task_states{idx_reach,2}(1:end, channel); 
                        vect   = [tmp_pre; tmp_post];  
                        firing_rate = [firing_rate, vect./bin_size];
                        
                    end
                end 
                zscored = (mean(firing_rate,2) - mean_baseline.by_targets{1,d}(target, ch_global))./std_baseline.by_targets{1,d}(target, ch_global);
                zscored(isnan(zscored) | isinf(zscored)) = 0;
                ZS_smooth = smoothdata(zscored, 1, 'gaussian', sigma_bins * 6);
                zscored_by_targets = [zscored_by_targets; ZS_smooth];   
            end 
            pca_matrix = [pca_matrix, zscored_by_targets]; 
        end
        pca_matrix_array = [pca_matrix_array, pca_matrix]; 
    end 

    condition = [condition; pca_matrix_array];
    condition_sep{d} = pca_matrix_array;
end 
% condition(:,14) = 0;
% condition(:,15) = 0; 
%% Filtraggio neuroni poco informativi
% meanFR = mean(condition, 1);          
% varFR = var(condition, 0, 1);
% 
% minMeanFR = 1;       % Hz: neuroni con FR medio < 1 Hz vengono scartati
% minVarFR  = 0.5;     % varianza minima (su firing rate in Hz)
% 
% valid_cols = (meanFR > minMeanFR) & (varFR > minVarFR);
% fprintf('\nNeuroni totali: %d, neuroni tenuti dopo filtro: %d\n', ...
%         numel(meanFR), sum(valid_cols));
% 
% condition = condition(:, valid_cols);
% for d = 1:numel(condition_sep)
%     condition_sep{d} = condition_sep{d}(:, valid_cols);
% end

%% Costruzione tensore per dPCA
% PCA matrix: [tempo × target × neuroni]
% X: [neurone, tempo, target, condizione]

n_bin = pre_bins + post_bins;
% n_bin = length(tmp_pre) + length(tmp_post); 
n_cond = numel(condition_sep);   
n_neurons = size(condition, 2);     

X = zeros(n_neurons, n_bin, n_targets, n_cond);
for c = 1:n_cond
    % Reshape: [tempo, target, neurone]
    cond_resh = reshape(condition_sep{c,1} , [n_bin, n_targets, n_neurons]);  
    
    % Permuta per ottenere [neurone, tempo, target]
    X(:,:,:,c) = permute(cond_resh, [3 1 2]);            
end

%% Vettore tempi (finestra PRE/POST)
%   -pre_bins*bin_size ... 0 ... (post_bins-1)*bin_size con t = 0 al Go

time = ((1:n_bin) - pre_bins - 1) * bin_size;   
timeEvents = 0;   % Go

%% Porta i dati nel formato dPCA
% firingRatesAverage: [N S D T]
% N = neuroni, S = target, D = condizione, T = tempo

N = size(X,1);    % 192 neuroni
T = size(X,2);    % time bins
S = size(X,3);    % 8 target
D = size(X,4);    % n condizioni

firingRatesAverage = permute(X, [1 3 4 2]);   % [N S D T]

%% Definizione marginalizzazioni 
% Parametri:
%   1 - target
%   2 - condizione
%   3 - tempo

combinedParams = { ...
   {1, [1 3]}, ...                 % 1: target-like (Goal)
   {2, [2 3]}, ...                 % 2: condition-like (Gaze status)
   {3}, ...                        % 3: time-only
   {[1 2], [1 2 3]} };             % 4: interazione target/condizione

   % {1, [1 3]}, ...                 % 1: target-like (Goal)
   % {2, [2 3]}, ...                 % 2: condition-like (Gaze status)

margNames   = {'Target', 'Gaze', 'Time', 'Tgt/Gaze Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

%% dPCA senza regularizzazione
Xmat = reshape(firingRatesAverage, size(firingRatesAverage,1), []); % N x (S*D*T)
varNeuron = var(Xmat, 0, 2);
keep = varNeuron > 1e-12;
firingRatesAverage = firingRatesAverage(keep,:,:,:);
N = sum(keep);
lambda = 1e-6;

nComponents = 30;  

[W, V, whichMarg] = dpca(firingRatesAverage, nComponents, ...
    'combinedParams', combinedParams, 'lambda', lambda);

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg, ...
    'time', time, ...
    'timeEvents', timeEvents, ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);

fprintf('Target-like:      %d\n', sum(whichMarg == 1));
fprintf('Condition-like:   %d\n', sum(whichMarg == 2));
fprintf('Time-only:        %d\n', sum(whichMarg == 3));
fprintf('Interaction T/C:  %d\n', sum(whichMarg == 4));

%% Componenti latenti Z (scores) calcolate da W
% Z: [nComponents x S x D x T]

Z = zeros(nComponents, S, D, T);

% Centro i dati sui neuroni
meanNeur = squeeze(mean(mean(mean(firingRatesAverage,4),3),2));   % [N x 1]
Xc = bsxfun(@minus, firingRatesAverage, reshape(meanNeur,[N 1 1 1]));  % [N S D T]

for k = 1:nComponents
    w = W(:,k);                              
    Zk = squeeze(sum(bsxfun(@times, Xc, reshape(w,[N 1 1 1])), 1));
    Z(k,:,:,:) = Zk;
end

%% Scelta delle componenti per la TDR 
%  2 assi "Goal" (target-like) + 1 asse "Gaze" (condition-like)

tgtMargID  = 1;   % Target-like 
condMargID = 2;   % Condition-like 

tgtIDs  = find(whichMarg == tgtMargID);
condIDs = find(whichMarg == condMargID);

% fallback basato sulla varianza tra target / condizioni
nComp = size(Z,1);
dirScore  = zeros(1,nComp);
condScore = zeros(1,nComp);

for k = 1:nComp
    Zk = squeeze(Z(k,:,:,:));    % [S x D x T]
    Zk = permute(Zk, [3 1 2]);   % [T x S x D]

    % varianza tra target (media su condizioni)
    mCond = mean(Zk, 3);         % [T x S]
    dirScore(k) = var(mCond(:));

    % varianza tra condizioni (media su target)
    mTgt = mean(Zk, 2);          % [T x 1 x D]
    mTgt = squeeze(mTgt);        % [T x D]
    condScore(k) = var(mTgt(:));
end

[~, dirOrder]  = sort(dirScore,  'descend');
[~, condOrder] = sort(condScore, 'descend');

if numel(tgtIDs) >= 2
    kGoal1 = tgtIDs(1);
    kGoal2 = tgtIDs(2);
else
    kGoal1 = dirOrder(1);
    kGoal2 = dirOrder(2);
end

if ~isempty(condIDs)
    kCondition = condIDs(1);
else
    kCondition = condOrder(1);
end

fprintf('Scelgo Goal dims = dPC %d, %d | Condition dim = dPC %d\n', ...
        kGoal1, kGoal2, kCondition);

%% Figure (1) - Plot 3D Target Direction Reduction (TDR)

Zg1 = squeeze(Z(kGoal1,:,:,:));     % [S x D x T]
Zg2 = squeeze(Z(kGoal2,:,:,:));     % [S x D x T]
Zp1 = squeeze(Z(kCondition,:,:,:)); % [S x D x T]

Zg1 = permute(Zg1, [3 1 2]);      % [T x S x D]
Zg2 = permute(Zg2, [3 1 2]);      % [T x S x D]
Zp1 = permute(Zp1, [3 1 2]);      % [T x S x D]

tMask = true(1, T);    

colors = [
    0.980, 0.525, 0.580;   % pink
    0.468, 0.751, 0.797;   % turchese/azzurro scuro
    0.311, 0.444, 0.506;   % blu molto scuro
    0.311, 0.444, 0.506;
    ];

legendLabels = cell(size(filename));
for i = 1:numel(filename)
    if contains(filename{i}, 'free-gaze', 'IgnoreCase', true)
        legendLabels{i} = 'free-gaze';
    elseif contains(filename{i}, 'motor', 'IgnoreCase', true)
        legendLabels{i} = 'gaze-on-center';
    elseif contains(filename{i}, 'controlled', 'IgnoreCase', true)
        legendLabels{i} = 'gaze-on-target';
    else
        legendLabels{i} = 'unknown';
    end
end

figure('Color','w'); 
hold on; grid on; box off;

w = 25; 
hCond = gobjects(D,1);   % handle per la legenda, uno per condizione

for cond = 1:D           % condizioni
    for tgt = 1:S        % target
        x = Zg1(tMask, tgt, cond);
        y = Zg2(tMask, tgt, cond);
        z = Zp1(tMask, tgt, cond);

        x = smoothdata(x, 'gaussian', w);
        y = smoothdata(y, 'gaussian', w);
        z = smoothdata(z, 'gaussian', w);

        h = plot3(x, y, z, 'LineWidth', 1.5, ...
                  'Color', colors(cond,:));  % curva

        if tgt == 1
            hCond(cond) = h;
        end

        plot3(x(1), y(1), z(1), 'o', ...
              'MarkerFaceColor', colors(cond,:), ...
              'MarkerEdgeColor', 'none', ...
              'MarkerSize', 6);

        plot3(x(end), y(end), z(end), '^', ...
              'MarkerFaceColor', colors(cond,:), ...
              'MarkerEdgeColor', 'none', ...
              'MarkerSize', 6);
    end
end

xlabel(sprintf('Goal Dim 1 (dPC %d)', kGoal1));
ylabel(sprintf('Goal Dim 2 (dPC %d)', kGoal2));
zlabel(sprintf('Gaze Dim 1 (dPC %d)', kCondition));

legend(hCond, legendLabels, 'Location', 'best');

view([-40 25]);
title('Subspace Goal/Gaze');

%% Figure (2) - Time vs Target vs Condition 

ciIDs = find(whichMarg == 3);           % condition-independent (tempo)
goalIDs = find(whichMarg == 1);           % gaze-independent (target)
gazeIDs = find(whichMarg == 2);           % main effect gaze
interIDs = find(whichMarg == 4);          % target×gaze

gazeDepIDs = [gazeIDs(:); interIDs(:)];      % tutto ciò che è gaze-dependent

compVar = explVar.componentVar;   % [1 x nComponents]

[~, idxCI] = max(compVar(ciIDs));    
kCI = ciIDs(idxCI);

[~, idxGoal] = max(compVar(goalIDs));
kGoal = goalIDs(idxGoal);

[~, idxGazeD] = max(compVar(gazeDepIDs));
kGazeDep = gazeDepIDs(idxGazeD);

Z_CI = permute(squeeze(Z(kCI    ,:,:,:)), [3 1 2]);  % [T x S x D]
Z_goal = permute(squeeze(Z(kGoal  ,:,:,:)), [3 1 2]);  % [T x S x D]
Z_gazeDp = permute(squeeze(Z(kGazeDep,:,:,:)), [3 1 2]); % [T x S x D]


colors_target = [
    0.839, 0.153, 0.157;  % rosso
    0.122, 0.467, 0.706;  % blu
    0.172, 0.627, 0.172;  % verde
    0.580, 0.404, 0.741;  % viola
    1.000, 0.498, 0.055;  % arancione
    0.737, 0.741, 0.133;  % giallo oliva
    0.549, 0.337, 0.294;  % marrone
    0.890, 0.466, 0.760;  % rosa
];

stylesGz  = {'-','--',':'};             % uno per ogni condizione di gaze (D = 3 max)

figure('Color','w'); 
w = 25; 

% Condition-independent (time)
subplot(1,3,1); hold on;
for d = 1:D          % condizioni di gaze
    for s = 1:S      % target
        z_ci = smoothdata(Z_CI(:,s,d), 'gaussian', w);
        
        plot(time, z_ci, ...
             'Color', colors_target(s,:), ...
             'LineStyle', stylesGz{d}, 'LineWidth', 1.2);
    end
end
title('Time-dependent');
xlabel('Time (s)'); ylabel('Projection onto dPC 1');
grid on;
box off;
axis tight;

% Gaze-independent (target)
subplot(1,3,2); hold on;
for d = 1:D
    for s = 1:S
        z_goal = smoothdata(Z_goal(:,s,d), 'gaussian', w);

        plot(time, z_goal, ...
             'Color', colors_target(s,:), ...
             'LineStyle', stylesGz{d}, 'LineWidth', 1.2);
    end
end
title('Target-dependent');
xlabel('Time (s)'); ylabel('Projection onto dPC 1');
grid on;
axis tight;

% Gaze-dependent
subplot(1,3,3); hold on;
for d = 1:D
    for s = 1:S
        z_gazeDp = smoothdata(Z_gazeDp(:,s,d), 'gaussian', w);

        plot(time, z_gazeDp, ...
             'Color', colors_target(s,:), ...
             'LineStyle', stylesGz{d}, 'LineWidth', 1.2);
    end
end
title('Gaze-dependent');
xlabel('Time (s)'); ylabel('Projection onto dPC 1');
grid on;
axis tight;

hold on;
hLeg = gobjects(D,1);
for d = 1:D
    hLeg(d) = plot(NaN, NaN, 'k', 'LineStyle', stylesGz{d}, 'LineWidth', 1.5);
end
legend(hLeg, legendLabels, 'Location', 'bestoutside');

%% Figure (3) - Istogramma
nShown = size(explVar.margVar,2);
margVar   = explVar.margVar(:, 1:nShown);   % [4 x nShown]
margTotal = sum(margVar, 2);                % [4 x 1]

% normalizzazione sulla varianza spiegata
pieVals   = 100 * margTotal / sum(margTotal);

% Seleziono Target (1), Condizione/Gaze (2), Time (3)
vals = pieVals(1:3);
newColors = [
    23 100 171;    % blu 
    187 20 25;     % rosso 
    150 150 150    % grigio 
] / 256;

labels = {'Target','Gaze','Time'};
figure('Color','w');
b = barh(vals, 'FaceColor', 'flat', 'EdgeColor','none');
b.CData = newColors;

clear set
set(gca, 'YTickLabel', labels);
xlabel('Explained variance (%)');

% Mostra i valori di percentuale accanto alle barre
for i = 1:length(vals)
    text(vals(i) + 0.5, i, sprintf('%.1f%%', vals(i)), ...
         'VerticalAlignment', 'middle');
end
box off;



%% Figure (4) - Fisso il target, vario la condizione
% Zcond = Zp1;           % [T x S x D]
% t = time(tMask);
% 
% if size(colors,1) < D
%     colorsCond = lines(D);
% else
%     colorsCond = colors(1:D,:);
% end
% 
% wSmooth = 25;
% 
% figure('Color','w'); 
% for tgt = 1:S
%     subplot(ceil(S/2), 2, tgt); hold on; box off;
%     for cond = 1:D
%         y = Zcond(tMask, tgt, cond);
%         y = smoothdata(y, 'gaussian', wSmooth);
%         plot(t, y, 'LineWidth', 1.5, ...
%              'Color', colorsCond(cond,:));
%     end
%     title(sprintf('Target %d', tgt), 'FontWeight','normal');
%     xlabel('Time (s)');
%     ylabel('dPC cond');
%     grid on;
%     axis tight;
% end
% legend(legendLabels, 'Position',[0.85 0.1 0.1 0.1]);
% % sgtitle('Gaze subspace');
% 
% 
% %% Figure (5) - Vario la condizione, fisso il target
% figure('Color','w'); 
% 
% for cond = 1:D
%     subplot(1, D, cond); hold on; box off;
%     for tgt = 1:S
%         y = Zcond(tMask, tgt, cond);
%         y = smoothdata(y, 'gaussian', wSmooth);
%         plot(t, y, 'LineWidth', 1.5, ...
%              'Color', colors_target(tgt,:));
%     end
%     title(sprintf('%s', legendLabels{cond}), 'FontWeight','normal');
%     xlabel('Time (s)');
%     ylabel('dPC cond');
%     grid on;
%     axis tight;
% end
% sgtitle('Gaze subspace');

