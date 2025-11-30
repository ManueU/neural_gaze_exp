% =========================================================
% DESCRIZIONE:
% Lo script esegue un’analisi dPCA (demixed PCA) su dati neurali
% provenienti dalla sola condizione "free-gaze", trattando i diversi
% set (es. 2,4,5,6) come condizioni sperimentali distinte.
%
% Obiettivi di questa versione:
%   - usare SOLO il dataset free-gaze;
%   - considerare ogni set come una "condizione" separata nella dPCA;
%   - verificare cosa succede alla TDR quando le condizioni dovrebbero
%     essere sostanzialmente identiche (attese traiettorie sovrapposte).
% =========================================================

% 1: medial arm 
% 2: lateral hand 

clearvars
close all
clc

sets = [1,2,4,5,6]; 
n_sets = numel(sets); 
n_arrays = 2;          
n_channels = 96;    
n_targets = 8; 
bin_size = 0.02;       
period_pre = 1.0;      
period_post = 0.8;     

% SOLO free-gaze
filename = '../00_Data_extraction/motor_BCI02.mat';

%% Caricamento dati e definizione finestra PRE/POST
fprintf('\nDataset (unico): %s\n', filename); 
load(filename);   

% stati PRE/POST per free-gaze
PRE  = "Pres12";
POST = "Reach";

% numero di bin PRE/POST
pre_bins  = max(1, round(period_pre/bin_size));
post_bins = max(1, round(period_post/bin_size));

%% Costruzione matrice condizione (media su trial, finestra PRE/POST)
% Ora OGNI SET è una condizione diversa: D = n_sets

condition    = []; 
condition_sep = cell(n_sets,1);   % una matrice per ogni set/condizione

for d = 1:n_sets
    setIdx = sets(d);
    fprintf('  -> Set come condizione: %d\n', setIdx);

    pca_matrix_array = []; 

    for array = 1:n_arrays
        pca_matrix = []; 

        % Trovo una volta gli indici degli stati PRE/POST per questo set/array
        idx_pres  = find(string(data(setIdx).Data(array).Interp(1).Task_states(:,1)) == PRE); 
        idx_reach = find(string(data(setIdx).Data(array).Interp(1).Task_states(:,1)) == POST); 
        if isempty(idx_pres) || isempty(idx_reach)
            error('Stati PRE/POST non trovati per set %d, array %d: controlla PRE/POST', ...
                  setIdx, array);
        end

        for channel = 1:n_channels
            firing_rate = []; 

            for target = 1:n_targets
                M_spikes = [];

                % selezione trial per target SOLO in questo set/condizione
                idx = find([data(setIdx).Data(array).Interp.Target_ID] == target); 

                for j = 1:length(idx)
                    % bin finestrati PRE/POST
                    tmp_pre_all  = data(setIdx).Data(array).Interp(idx(j)).Task_states{idx_pres, 2}(:, channel); 
                    tmp_post_all = data(setIdx).Data(array).Interp(idx(j)).Task_states{idx_reach,2}(:, channel); 

                    if size(tmp_pre_all,1) < pre_bins
                        error('Trial troppo corto nella finestra PRE (set %d, array %d).', setIdx, array);
                    end
                    if size(tmp_post_all,1) < post_bins
                        error('Trial troppo corto nella finestra POST (set %d, array %d).', setIdx, array);
                    end

                    tmp_pre  = tmp_pre_all(end-pre_bins+1:end);
                    tmp_post = tmp_post_all(1:post_bins);

                    matrix   = [tmp_pre; tmp_post];  
                    M_spikes = [M_spikes, matrix];   
                end

                % media sui trial per target, poi FR in Hz
                M_spikes_mean = mean(M_spikes, 2);     
                firing_rate   = [firing_rate; M_spikes_mean ./ bin_size]; 
            end 

            % concatena i diversi canali
            pca_matrix = [pca_matrix, firing_rate]; 
        end

        % concatena i due array
        pca_matrix_array = [pca_matrix_array, pca_matrix]; 
    end 

    % ogni set è una "condizione" separata
    condition = [condition; pca_matrix_array];
    condition_sep{d} = pca_matrix_array;
end 

%% Filtraggio neuroni poco informativi
meanFR = mean(condition, 1);          
varFR  = var(condition, 0, 1);

minMeanFR = 1;       % Hz: neuroni con FR medio < 1 Hz vengono scartati
minVarFR  = 0.5;     % varianza minima (su firing rate in Hz)

valid_cols = (meanFR > minMeanFR) & (varFR > minVarFR);
fprintf('\nNeuroni totali: %d, neuroni tenuti dopo filtro: %d\n', ...
        numel(meanFR), sum(valid_cols));

condition = condition(:, valid_cols);
for d = 1:numel(condition_sep)
    condition_sep{d} = condition_sep{d}(:, valid_cols);
end

%% Costruzione tensore per dPCA
% X: [neurone, tempo, target, condizione]

n_bin     = pre_bins + post_bins;
n_cond    = numel(condition_sep);      % ora = n_sets
n_neurons = size(condition, 2);     

X = zeros(n_neurons, n_bin, n_targets, n_cond);
for d = 1:n_cond
    % Reshape: [tempo, target, neurone]
    cond_resh = reshape(condition_sep{d}, [n_bin, n_targets, n_neurons]);  
    
    % Permuta per ottenere [neurone, tempo, target]
    X(:,:,:,d) = permute(cond_resh, [3 1 2]);            
end

%% Vettore tempi (finestra PRE/POST)
%   -pre_bins*bin_size ... 0 ... (post_bins-1)*bin_size con t = 0 al Go

time      = ((1:n_bin) - pre_bins - 1) * bin_size;   
timeEvents = 0;   % Go

%% Porta i dati nel formato dPCA
% firingRatesAverage: [N S D T]

N = size(X,1);    
T = size(X,2);    
S = size(X,3);    
D = size(X,4);    % ora D = n_sets (numero di set usati)

firingRatesAverage = permute(X, [1 3 4 2]);   % [N S D T]

%% Definizione marginalizzazioni 
% Parametri:
%   1 - target
%   2 - condizione (qui: SET)
%   3 - tempo

combinedParams = { ...
   {1, [1 3]}, ...                 % 1: target-like (Goal)
   {2, [2 3]}, ...                 % 2: condition-like (Set)
   {3}, ...                        % 3: time-only
   {[1 2], [1 2 3]} };             % 4: interazione target/condizione

margNames   = {'Target', 'Set', 'Time', 'Tgt/Set Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

%% dPCA senza regularizzazione
nComponents = 30;  

[W, V, whichMarg] = dpca(firingRatesAverage, nComponents, ...
    'combinedParams', combinedParams);

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

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
%  2 assi "Goal" (target-like) + 1 asse "Condition-like" (qui: set)

tgtMargID  = 1;   % Target-like 
condMargID = 2;   % Condition-like 

tgtIDs  = find(whichMarg == tgtMargID);
condIDs = find(whichMarg == condMargID);

% fallback basato sulla varianza tra target / condizioni
nComp     = size(Z,1);
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

% colori: uno per ogni set/condizione
colors = lines(D);

% etichette leggenda: un nome per ogni set
legendLabels = arrayfun(@(s) sprintf('set %d', s), sets, 'UniformOutput', false);

figure('Color','w'); 
hold on; grid on; box off;

w = 25; 
hCond = gobjects(D,1);   % handle per la legenda, uno per condizione

for cond = 1:D           % condizioni (qui: set)
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
zlabel(sprintf('Condition Dim 1 (dPC %d)', kCondition));

legend(hCond, legendLabels, 'Location', 'best');

view([-40 25]);
title('Subspace Goal / Set (free-gaze only)');

%% Figure (2) - Time vs Target vs Condition 

ciIDs    = find(whichMarg == 3);  % time
goalIDs  = find(whichMarg == 1);  % target
gazeIDs  = find(whichMarg == 2);  % qui: set
interIDs = find(whichMarg == 4);  % target×set

gazeDepIDs = [gazeIDs(:); interIDs(:)];      % tutto ciò che è condition-dependent

compVar = explVar.componentVar;   % [1 x nComponents]

[~, idxCI]    = max(compVar(ciIDs));    
kCI           = ciIDs(idxCI);

[~, idxGoal]  = max(compVar(goalIDs));
kGoal         = goalIDs(idxGoal);

[~, idxGazeD] = max(compVar(gazeDepIDs));
kGazeDep      = gazeDepIDs(idxGazeD);

Z_CI     = permute(squeeze(Z(kCI    ,:,:,:)), [3 1 2]);  % [T x S x D]
Z_goal   = permute(squeeze(Z(kGoal  ,:,:,:)), [3 1 2]);  % [T x S x D]
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

% uno stile per condizione/set (qui possono essere di più, ne mettiamo 4)
stylesGz  = {'-','--',':','-.'};

figure('Color','w'); 
w = 25; 

% Condition-independent (time)
subplot(1,3,1); hold on;
for d = 1:D          % condizioni (set)
    for s = 1:S      % target
        z_ci = smoothdata(Z_CI(:,s,d), 'gaussian', w);
        
        plot(time, z_ci, ...
             'Color', colors_target(s,:), ...
             'LineStyle', stylesGz{min(d,numel(stylesGz))}, ...
             'LineWidth', 1.2);
    end
end
title('Time-dependent');
xlabel('Time (s)'); ylabel('Projection onto dPC 1');
grid on;
box off;
axis tight;

% Target-dependent
subplot(1,3,2); hold on;
for d = 1:D
    for s = 1:S
        z_goal_  = smoothdata(Z_goal(:,s,d), 'gaussian', w);

        plot(time, z_goal_, ...
             'Color', colors_target(s,:), ...
             'LineStyle', stylesGz{min(d,numel(stylesGz))}, ...
             'LineWidth', 1.2);
    end
end
title('Target-dependent');
xlabel('Time (s)'); ylabel('Projection onto dPC 1');
grid on;
axis tight;

% Condition-dependent (Set)
subplot(1,3,3); hold on;
for d = 1:D
    for s = 1:S
        z_gazeDp_ = smoothdata(Z_gazeDp(:,s,d), 'gaussian', w);

        plot(time, z_gazeDp_, ...
             'Color', colors_target(s,:), ...
             'LineStyle', stylesGz{min(d,numel(stylesGz))}, ...
             'LineWidth', 1.2);
    end
end
title('Set-dependent');
xlabel('Time (s)'); ylabel('Projection onto dPC 1');
grid on;
axis tight;

hold on;
hLeg = gobjects(D,1);
for d = 1:D
    hLeg(d) = plot(NaN, NaN, 'k', ...
                   'LineStyle', stylesGz{min(d,numel(stylesGz))}, ...
                   'LineWidth', 1.5);
end
legend(hLeg, legendLabels, 'Location', 'bestoutside');

%% Figure (3) - Istogramma
nShown    = size(explVar.margVar,2);
margVar   = explVar.margVar(:, 1:nShown);   % [4 x nShown]
margTotal = sum(margVar, 2);                % [4 x 1]

% normalizzazione sulla varianza spiegata
pieVals = 100 * margTotal / sum(margTotal);

% Seleziono Target (1), Condition/Set (2), Time (3)
vals = pieVals(1:3);
newColors = [
    23 100 171;    % blu 
    187 20 25;     % rosso 
    150 150 150    % grigio 
] / 256;

labels = {'Target','Set','Time'};
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

%% Figure (4) - Fisso il target, vario la condizione (set)
Zcond = Zp1;           % [T x S x D]
t     = time(tMask);

if size(colors,1) < D
    colorsCond = lines(D);
else
    colorsCond = colors(1:D,:);
end

wSmooth = 25;

figure('Color','w'); 
for tgt = 1:S
    subplot(ceil(S/2), 2, tgt); hold on; box off;
    for cond = 1:D
        y = Zcond(tMask, tgt, cond);
        y = smoothdata(y, 'gaussian', wSmooth);
        plot(t, y, 'LineWidth', 1.5, ...
             'Color', colorsCond(cond,:));
    end
    title(sprintf('Target %d', tgt), 'FontWeight','normal');
    xlabel('Time (s)');
    ylabel('dPC cond');
    grid on;
    axis tight;
end
legend(legendLabels, 'Position',[0.85 0.1 0.1 0.1]);

%% Figure (5) - Vario la condizione (set), fisso il target
figure('Color','w'); 

for cond = 1:D
    subplot(1, D, cond); hold on; box off;
    for tgt = 1:S
        y = Zcond(tMask, tgt, cond);
        y = smoothdata(y, 'gaussian', wSmooth);
        plot(t, y, 'LineWidth', 1.5, ...
             'Color', colors_target(tgt,:));
    end
    title(sprintf('%s', legendLabels{cond}), 'FontWeight','normal');
    xlabel('Time (s)');
    ylabel('dPC target');
    grid on;
    axis tight;
end
