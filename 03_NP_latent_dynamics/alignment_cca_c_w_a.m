% =========================================================
% DESCRIZIONE:
% Implementa un’analisi ispirata a Safaie et al. (2023, Nature) per
% quantificare quanto le dinamiche neurali siano preservate:
% - tra due condizioni sperimentali diverse (ACROSS conditions);
% - all’interno della stessa condizione usando uno split trial-by-trial
%   (WITHIN, upper bound);
% - tra finestre di controllo non allineate al comportamento (CONTROL,
%   lower bound).
%
% L’analisi si svolge nello spazio latente comune ottenuto tramite PCA
% applicata a tutte le condizioni, e confronta le traiettorie neurali
% per ciascun target tramite Canonical Correlation Analysis (CCA).
%
% La procedura prevede:
% 1) estrazione, per ciascun dataset, array, canale e target, delle
%    risposte neurali in una finestra temporale fissa (reach), 
%    combinando periodo pre- e post-evento;
% 2) costruzione di:
%       - matrici “behavioural” (condition_all) per PCA;
%       - split A/B dei trial (condition_A_all, condition_B_all) per la
%         CCA WITHIN;
%       - finestre di controllo random (condition_ctrl_all) per la CCA CONTROL;
% 3) media delle risposte per target e conversione in firing rate
%    (divisione per bin_size);
% 4) concatenazione tra array e tra condizioni, con rimozione delle
%    colonne a varianza nulla;
% 5) z-score globale e PCA sullo spazio comune (coeff, score, explained);
% 6) visualizzazione delle traiettorie nelle prime 3 PC per tutte le
%    condizioni e target (Figura 1);
% 7) CCA ACROSS: confronto tra le due condizioni principali usando le
%    prime PC (pc_idx_cca) come input (Figura 2);
% 8) CCA WITHIN: confronto tra split A/B della condizione 1 (upper bound);
% 9) CCA CONTROL: confronto tra finestre di controllo di condizione 1 e 2
%    con shuffle temporale delle righe (lower bound);
% 10) confronto finale dei valori di correlazione canonica ACROSS,
%     WITHIN e CONTROL lungo le dimensioni canoniche (Figura 3).
%
% Prima dell’esecuzione, assicurarsi di impostare correttamente:
%   filename    : cell array con i percorsi ai file .mat delle condizioni
%   sets        : indici dei set/trial block da includere
%   n_arrays    : numero di array registrati
%   n_channels  : numero di canali per array
%   n_targets   : numero di target comportamentali
%   bin_size    : dimensione del bin temporale [s]
%   period_pre  : finestra pre-evento rispetto allo stato PRE [s]
%   period_post : finestra post-evento rispetto allo stato POST [s]
%   PRE, POST   : etichette degli stati comportamentali per ciascun dataset
%                 (derivate da Task_states e impostate nello script)
% 
% REF: Safaie, M., Chang, J. C., Park, J., Miller, L. E., Dudman, J. T., Perich, M. G., 
% & Gallego, J. A. (2023). Preserved neural dynamics across animals performing similar 
% behaviour. Nature, 623(7988), 765-771.
% =========================================================


%% CCA across / within / control

clearvars
% close all
clc

%% Parametri di base
sets       = [1,2,4,5,6]; 
n_sets     = numel(sets); 
n_arrays   = 2;      
n_channels = 96;    
n_targets  = 8; 

bin_size    = 0.02;   % [s]
period_pre  = 0.1;    % [s] finestra prima dell'evento
period_post = 0.5;    % [s] finestra dopo l'evento

pre_bins  = max(1, round(period_pre/bin_size));
post_bins = max(1, round(period_post/bin_size));
nbin      = pre_bins + post_bins;

filename = { ...
    '../00_Data_extraction/controlled_BCI02.mat', ...
    '../00_Data_extraction/motor_BCI02.mat'};

nCond = numel(filename);

%% Matrici per PCA, within e control
condition_all      = [];                    % tutte le condizioni, finestra "comportamentale"
condition_A_all    = [];                    % split A per tutte le condizioni
condition_B_all    = [];                    % split B per tutte le condizioni
condition_ctrl_all = [];                    % finestre di controllo per tutte le condizioni

condition_sep_all  = cell(nCond,1);         % per condizione (behavioural)
condition_sep_A    = cell(nCond,1);         % per condizione (within A)
condition_sep_B    = cell(nCond,1);         % per condizione (within B)
condition_sep_ctrl = cell(nCond,1);         % per condizione (control windows)


%% Costruzione delle matrici
for d = 1:nCond
    
    fprintf('\nDataset: %s\n', filename{d}); 
    load(filename{d});  
    
    [~, baseName, ext] = fileparts(filename{d});
    ds_name = [baseName ext];

    % Stati PRE/POST a seconda del dataset
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
        
    idx_pres  = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == PRE); 
    idx_reach = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == POST); 
    if isempty(idx_pres) || isempty(idx_reach)
        error('Stati PRE/POST non trovati: controlla PRE/POST per il dataset %s', ds_name);
    end

    % Matrici per il dataset d-esimo
    pca_matrix_array_all  = []; 
    pca_matrix_array_A    = [];
    pca_matrix_array_B    = [];
    pca_matrix_array_ctrl = [];
    
    for array = 1:n_arrays
        
        pca_matrix_all  = []; 
        pca_matrix_A    = [];
        pca_matrix_B    = [];
        pca_matrix_ctrl = [];
    
        for channel = 1:n_channels
            
            firing_rate_all  = []; 
            firing_rate_A    = [];
            firing_rate_B    = [];
            firing_rate_ctrl = [];
            
            for target = 1:n_targets
                
                M_spikes_all  = [];
                M_spikes_A    = [];
                M_spikes_B    = [];
                M_spikes_ctrl = [];
    
                trial_counter = 0;
    
                for set_ = 1:n_sets
                    set = sets(set_);
                    idx_trials = find([data(set).Data(array).Resampled.Target_ID] == target); 
                    
                    for j = 1:length(idx_trials)
                        trial_counter = trial_counter + 1;
                        tr = idx_trials(j);

                        % ---------- finestra COMPORTAMENTALE (pre+post) ----------
                        tmp_pre  = data(set).Data(array).Resampled(tr).Task_states{idx_pres,  2}(end-pre_bins+1:end, channel); 
                        tmp_post = data(set).Data(array).Resampled(tr).Task_states{idx_reach, 2}(1:post_bins,        channel); 
                        matrix   = [tmp_pre; tmp_post];      % (nbin x 1)
                        
                        M_spikes_all = [M_spikes_all, matrix];

                        % split A/B (within)
                        if mod(trial_counter,2) == 1
                            M_spikes_A = [M_spikes_A, matrix];
                        else
                            M_spikes_B = [M_spikes_B, matrix];
                        end

                        % ---------- finestra di CONTROLLO  ----------
                        % concateno tutti gli stati del trial (trial+intertrial)
                        trial_full = [];
                        nStates = size(data(set).Data(array).Resampled(tr).Task_states,1);
                        for s = 1:nStates
                            trial_full = [trial_full; ...
                                data(set).Data(array).Resampled(tr).Task_states{s,2}(:, channel)];
                        end

                        L = size(trial_full,1);
                        if L >= nbin
                            start_idx = randi([1, L - nbin + 1]);
                            ctrl_seg  = trial_full(start_idx:start_idx+nbin-1);
                        else
                            % fallback: se il trial è troppo corto, uso la finestra comportamentale
                            ctrl_seg = matrix;
                        end
                        M_spikes_ctrl = [M_spikes_ctrl, ctrl_seg];
                    end
                end 

                % ---------- medie per target ----------
                M_spikes_mean_all = mean(M_spikes_all, 2);

                if ~isempty(M_spikes_A)
                    M_spikes_mean_A = mean(M_spikes_A, 2);
                else
                    M_spikes_mean_A = M_spikes_mean_all;
                end

                if ~isempty(M_spikes_B)
                    M_spikes_mean_B = mean(M_spikes_B, 2);
                else
                    M_spikes_mean_B = M_spikes_mean_all;
                end

                if ~isempty(M_spikes_ctrl)
                    M_spikes_mean_ctrl = mean(M_spikes_ctrl, 2);
                else
                    M_spikes_mean_ctrl = M_spikes_mean_all;
                end
    
                firing_rate_all  = [firing_rate_all;  M_spikes_mean_all  ./ bin_size];
                firing_rate_A    = [firing_rate_A;    M_spikes_mean_A    ./ bin_size];
                firing_rate_B    = [firing_rate_B;    M_spikes_mean_B    ./ bin_size];
                firing_rate_ctrl = [firing_rate_ctrl; M_spikes_mean_ctrl ./ bin_size];
    
            end % target

            % un canale = una colonna
            pca_matrix_all  = [pca_matrix_all,  firing_rate_all]; 
            pca_matrix_A    = [pca_matrix_A,    firing_rate_A];
            pca_matrix_B    = [pca_matrix_B,    firing_rate_B];
            pca_matrix_ctrl = [pca_matrix_ctrl, firing_rate_ctrl];
    
        end % channel

        % concate dei due array
        pca_matrix_array_all  = [pca_matrix_array_all,  pca_matrix_all]; 
        pca_matrix_array_A    = [pca_matrix_array_A,    pca_matrix_A];
        pca_matrix_array_B    = [pca_matrix_array_B,    pca_matrix_B];
        pca_matrix_array_ctrl = [pca_matrix_array_ctrl, pca_matrix_ctrl];
    
    end % array

    % accodo tra condizioni
    condition_all      = [condition_all;      pca_matrix_array_all];
    condition_A_all    = [condition_A_all;    pca_matrix_array_A];
    condition_B_all    = [condition_B_all;    pca_matrix_array_B];
    condition_ctrl_all = [condition_ctrl_all; pca_matrix_array_ctrl];

    condition_sep_all{d}  = pca_matrix_array_all;
    condition_sep_A{d}    = pca_matrix_array_A;
    condition_sep_B{d}    = pca_matrix_array_B;
    condition_sep_ctrl{d} = pca_matrix_array_ctrl;
end 


%% Rimozione colonne a varianza nulla
col_std   = std(condition_all, 0, 1);
valid_cols = col_std > 0;

condition_all      = condition_all(:,      valid_cols);
condition_A_all    = condition_A_all(:,    valid_cols);
condition_B_all    = condition_B_all(:,    valid_cols);
condition_ctrl_all = condition_ctrl_all(:, valid_cols);

for d = 1:nCond
    condition_sep_all{d}  = condition_sep_all{d}(:,  valid_cols);
    condition_sep_A{d}    = condition_sep_A{d}(:,    valid_cols);
    condition_sep_B{d}    = condition_sep_B{d}(:,    valid_cols);
    condition_sep_ctrl{d} = condition_sep_ctrl{d}(:, valid_cols);
end


%% PCA sullo spazio comune
[Xz, muZ, sigmaZ] = zscore(condition_all, 0, 1); 
[coeff, score, latent, tsquared, explained] = pca(Xz, 'Algorithm','svd');


%% Figura (1) - PCA (trajectories PC1-3)
T      = repmat(n_targets * nbin, nCond, 1); 
offset = [0; cumsum(T(1:end-1))];

figure('Color','White'); hold on
styles = {'-','--',':','-.'};

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

w = 11; 

for c = 1:nCond
    for t = 1:n_targets
        rows = offset(c) + (t-1)*nbin + (1:nbin);
        traj = score(rows,1:3);
        t1 = smoothdata(traj(:,1), 'gaussian', w);
        t2 = smoothdata(traj(:,2), 'gaussian', w);
        t3 = smoothdata(traj(:,3), 'gaussian', w);

        plot3(t1, t2, t3, 'Color', colors(t,:), 'LineWidth', 2, ...
              'LineStyle', styles{mod(c-1,numel(styles))+1});
    
        plot3(t1(1), t2(1), t3(1), 'o', ...
              'MarkerSize', 8, 'MarkerFaceColor', colors(t,:), ...
              'MarkerEdgeColor', colors(t,:), 'HandleVisibility','off');
    
        plot3(t1(end), t2(end), t3(end), '^', ...
              'MarkerSize', 8, 'MarkerFaceColor', colors(t,:), ...
              'MarkerEdgeColor', colors(t,:), 'HandleVisibility','off');
    end
end 

hStart = plot3(nan,nan,nan, 'o', 'Color','k', 'MarkerFaceColor','k');
hEnd   = plot3(nan,nan,nan, '^', 'Color','k', 'MarkerFaceColor','k');

if nCond == 1
    legend([hStart hEnd], {'Start','End'}, 'Location','northeastoutside');
else
    legEntries = cell(1, nCond + 2);  
    styles = {'-','--',':','-.'};
    hStyles = gobjects(1, nCond);

    for c = 1:nCond
        [~, baseName, ~] = fileparts(filename{c});
        [namePart, ~] = strtok(baseName, '_');
        legEntries{c} = namePart;
        ls = styles{mod(c-1, numel(styles))+1};
        hStyles(c) = plot3(nan,nan,nan, ls, 'Color','k', 'LineWidth', 2);
    end

    legEntries{nCond+1} = 'Start';
    legEntries{nCond+2} = 'End';
    legend([hStyles hStart hEnd], legEntries, 'Location','northeastoutside');
end

xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
grid on; axis equal; 
title('Traiettorie PCA');


%% CCA ACROSS condizioni
Xz_cond1 = (condition_sep_all{1} - muZ) ./ sigmaZ;
Xz_cond2 = (condition_sep_all{2} - muZ) ./ sigmaZ;

scores1 = Xz_cond1 * coeff;   
scores2 = Xz_cond2 * coeff;

pc_idx_cca = 1:10;

Xcca_across = zscore(scores1(:, pc_idx_cca));
Ycca_across = zscore(scores2(:, pc_idx_cca));

[A_across, B_across, r_across, U_across, V_across, stats_across] = canoncorr(Xcca_across, Ycca_across); 
disp('Canonical correlations (CCA ACROSS condizioni):');
disp(r_across);


%% Figure (2) - Traiettorie nello spazio canonico (across)
figure('Color','w'); hold on;
axis equal; grid on;
xlabel('Can1'); ylabel('Can2'); zlabel('Can3');
title('Traiettorie nello spazio canonico (ACROSS)');

for tg = 1:n_targets
    idx = (tg-1)*nbin + (1:nbin);

    Ut = U_across(idx, 1:3);   % cond1
    Vt = V_across(idx, 1:3);   % cond2

    u1 = smoothdata(Ut(:,1), 'gaussian', w);
    u2 = smoothdata(Ut(:,2), 'gaussian', w);
    u3 = smoothdata(Ut(:,3), 'gaussian', w);

    v1 = smoothdata(Vt(:,1), 'gaussian', w);
    v2 = smoothdata(Vt(:,2), 'gaussian', w);
    v3 = smoothdata(Vt(:,3), 'gaussian', w);

    plot3(u1, u2, u3, '-', ...
        'Color', colors(tg,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Target %d - cond1', tg));

    plot3(v1, v2, v3, '--', ...
        'Color', colors(tg,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Target %d - cond2', tg));
end

h1 = plot3(nan,nan,nan,'-k','LineWidth',2,'DisplayName','Condizione 1');
h2 = plot3(nan,nan,nan,'--k','LineWidth',2,'DisplayName','Condizione 2');
legend([h1 h2],'Location','northeastoutside');


%% CCA WITHIN (upper bound) sulla condizione 1
Xz_cond1_A = (condition_sep_A{1} - muZ) ./ sigmaZ;
Xz_cond1_B = (condition_sep_B{1} - muZ) ./ sigmaZ;

scores1_A = Xz_cond1_A * coeff;
scores1_B = Xz_cond1_B * coeff;

X1 = zscore(scores1_A(:, pc_idx_cca));
X2 = zscore(scores1_B(:, pc_idx_cca));

[A_w, B_w, r_within, U_within, V_within] = canoncorr(X1, X2); 
disp('Canonical correlations WITHIN condizione 1 (upper bound):');
disp(r_within);


%% CCA CONTROL (lower bound)
Xz_ctrl1 = (condition_sep_ctrl{1} - muZ) ./ sigmaZ;
Xz_ctrl2 = (condition_sep_ctrl{2} - muZ) ./ sigmaZ;

scores_ctrl1 = Xz_ctrl1 * coeff;
scores_ctrl2 = Xz_ctrl2 * coeff;

X_ctrl = zscore(scores_ctrl1(:, pc_idx_cca));
Y_ctrl = zscore(scores_ctrl2(:, pc_idx_cca));

% shuffle delle righe per rompere la corrispondenza target×tempo
perm_rows  = randperm(size(Y_ctrl,1));
Y_ctrl_shuf = Y_ctrl(perm_rows, :);

[~, ~, r_control] = canoncorr(X_ctrl, Y_ctrl_shuf);
disp('Canonical correlations CONTROL (lower bound):');
disp(r_control);


%% Figure (3) - Figura finale: Across vs Within vs Control
figure('Color','w'); 
hold on;

plot(r_across,  '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
    'Color', [0.2 0.6 0.2], 'DisplayName', 'Across conditions');

plot(r_within,  '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
    'Color', [0.2 0.2 0.8], 'DisplayName', 'Within (upper bound)');

plot(r_control, '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
    'Color', [0.5 0.5 0.5], 'DisplayName', 'Control (lower bound)');

ylim([0 1])
xlabel('Canonical dimension');
ylabel('Canonical correlation');
legend('Location', 'northeast');
grid on; box off;
