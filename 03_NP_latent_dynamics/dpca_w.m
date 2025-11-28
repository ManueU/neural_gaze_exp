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
period_pre = 0.05;      
period_post = 0.2;        

filename = { ...
    '../00_Data_extraction/controlled_BCI02.mat', ...
    '../00_Data_extraction/motor_BCI02.mat'};   


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
    pre_bins  = max(1, round(period_pre/bin_size));
    post_bins = max(1, round(period_post/bin_size));
    
    pca_matrix_array = []; 
    for array = 1:n_arrays
        pca_matrix = []; 
    
        for channel = 1:n_channels
            firing_rate = []; 
            
            for target = 1:n_targets
                M_spikes = [];
    
                for set_ = 1:n_sets
                    set = sets(set_);
                    idx = find([data(set).Data(array).Interp.Target_ID] == target); 
    
                    for j = 1:length(idx)
                        % bin finestrati PRE/POST
                        tmp_pre  = data(set).Data(array).Interp(idx(j)).Task_states{idx_pres, 2}(end-pre_bins+1:end, channel); 
                        tmp_post = data(set).Data(array).Interp(idx(j)).Task_states{idx_reach,2}(1:post_bins, channel); 
                        matrix   = [tmp_pre; tmp_post];  
                        
                        M_spikes = [M_spikes, matrix];   
                    end
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

    condition = [condition; pca_matrix_array];
    condition_sep{d} = pca_matrix_array;
end 

%% Filtraggio neuroni poco informativi
meanFR = mean(condition, 1);          
varFR = var(condition, 0, 1);

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

n_bin = pre_bins + post_bins;
n_cond = numel(condition_sep);   
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

time = ((1:n_bin) - pre_bins - 1) * bin_size;   
timeEvents = 0;   % Go

%% Porta i dati nel formato dPCA
% firingRatesAverage: [N S D T]
% N = neuroni, S = target, D = condizione, T = tempo

N = size(X,1);    % 192 neuroni
T = size(X,2);    % time bins
S = size(X,3);    % 8 target
D = size(X,4);    % 2 condizioni

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

margNames   = {'Target', 'Gaze', 'Time', 'Tgt/Gaze Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

%% dPCA senza regularizzazione
nComponents = 30;  

[W, V, whichMarg] = dpca(firingRatesAverage, nComponents, ...
    'combinedParams', combinedParams);

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

% centro i dati sui neuroni (come fa dPCA internamente)
meanNeur = squeeze(mean(mean(mean(firingRatesAverage,4),3),2));   % [N x 1]
Xc = bsxfun(@minus, firingRatesAverage, reshape(meanNeur,[N 1 1 1]));  % [N S D T]

for k = 1:nComponents
    w = W(:,k);                              
    % proiezione: somma sui neuroni pesata da w
    % risultato: [S x D x T]
    Zk = squeeze(sum(bsxfun(@times, Xc, reshape(w,[N 1 1 1])), 1));
    Z(k,:,:,:) = Zk;                         % [1 S D T]
end

%% Scelta delle componenti per la TDR 
%  2 assi "Goal" (target-like) + 1 asse "Posture" (condition-like)

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

%% Plot 3D stile TDR: Free-gaze vs Controlled
% Prendiamo le traiettorie dai punteggi Z
Zg1 = squeeze(Z(kGoal1,:,:,:));   % [S x D x T]
Zg2 = squeeze(Z(kGoal2,:,:,:));   % [S x D x T]
Zp1 = squeeze(Z(kCondition,:,:,:)); % [S x D x T]

% portiamo a [T x S x D]
Zg1 = permute(Zg1, [3 1 2]);      % [T x S x D]
Zg2 = permute(Zg2, [3 1 2]);      % [T x S x D]
Zp1 = permute(Zp1, [3 1 2]);      % [T x S x D]

% Il segnale è finestrato PRE/POST: uso TUTTI i bin
tMask = true(1, T);    
colors = lines(D);     % D=2: 1=free-gaze, 2=controlled

figure('Color','w'); 
hold on; 
grid on; 
box off;

for cond = 1:D           % 1 = free-gaze, 2 = controlled
    for tgt = 1:S        % 8 target
        x = Zg1(tMask, tgt, cond);
        y = Zg2(tMask, tgt, cond);
        z = Zp1(tMask, tgt, cond);

        % traiettoria 3D nel subspazio Goal/Gaze
        plot3(x, y, z, 'LineWidth', 1.5, ...
              'Color', colors(cond,:));

        % punto al Go (tempo 0)
        [~,i0] = min(abs(time));   % indice Go
        i0 = max(1, min(sum(tMask), i0));  % sicurezza
        plot3(x(i0), y(i0), z(i0), 'o', ...
              'MarkerFaceColor', colors(cond,:), ...
              'MarkerEdgeColor', 'none', ...
              'MarkerSize', 6);

        % punto alla fine della finestra (ultimo bin usato)
        plot3(x(end), y(end), z(end), 'o', ...
              'MarkerFaceColor', 'k', ...
              'MarkerEdgeColor', 'none', ...
              'MarkerSize', 4);
    end
end

xlabel(sprintf('Goal Dim 1 (dPC %d)', kGoal1));
ylabel(sprintf('Goal Dim 2 (dPC %d)', kGoal2));
zlabel(sprintf('Gaze Dim 1 (dPC %d)', kCondition));

% legenda "Gaze 1 vs Gaze 2"
hFG = plot3(nan,nan,nan,'-','Color',colors(1,:),'LineWidth',1.5);
hCT = plot3(nan,nan,nan,'-','Color',colors(2,:),'LineWidth',1.5);
legend([hFG hCT], {'Free-gaze','Gaze-constrained'}, 'Location','northeastoutside');

view([-40 25]);
title('Subspazio Goal/Gaze');
