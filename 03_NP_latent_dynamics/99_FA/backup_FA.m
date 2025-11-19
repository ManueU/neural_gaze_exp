clear all
close all
clc 

period_pre = 0.1; 
period_reach = 0.5;
n_trials = 32; % trials per set 
bin_size = 0.02; 
n_targets = 8; 

 
mat_files = { ...
    'motor_BCI02.mat' ... 
    'free-gaze_BCI02.mat'
};



%% FA on single condition
% sigma = 0.2; 
% window = round(5 * sigma / bin_size);  % perché σ ≈ window/5
% window = max(3, 2*floor(window/2)+1);
for d = 1:numel(mat_files) 
    disp(mat_files(d)); 
    ds_name = mat_files{d};
    load(ds_name);

    % Matrix construction
    if strcmp(ds_name, 'controlled_BCI02.mat')
        PRES = "Pres12"; 
        REACH = "Gaze"; 
    end 

    if strcmp(ds_name, 'gaze_BCI02.mat')
        PRES = "Pres12"; 
        REACH = "Gaze";
    else 
        PRES = "Pres12";
        REACH = "Reach";
    end 

    idx_pres = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == PRES); 
    idx_reach = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == REACH); 
    
    start_pres = size(data(1).Data(2).Resampled(1).Task_states{idx_pres,2},1) - round(period_pre/bin_size) + 1;
    end_reach  = round(period_reach/bin_size);
    
    pca_matrix_array = [];
    for array = 1:2
        pca_matrix = [];

        for channel = 1:96
            firing_rate = [];

            for target = 1:n_targets
                M_spikes_tmp = [];
        
                for set = 1:6
                    idx = find([data(set).Data(array).Resampled.Target_ID] == target); 
    
                    for j = 1:length(idx)
                        tmp_pres = data(set).Data(array).Resampled(idx(j)).Task_states{idx_pres, 2}(start_pres:end, channel); 
                        tmp_reach = data(set).Data(array).Resampled(idx(j)).Task_states{idx_reach,2}(1:end_reach, channel); 
                        matrix = [tmp_pres; tmp_reach]; 
                        
                        % Square root transformation
                        % matrix = sqrt(matrix); 
                        % Gaussian smoothing    
                        % m_trial = smoothdata(m_trial ./ bin_size, 1, 'gaussian', window); 
                        
                        M_spikes = [M_spikes_tmp, matrix];
                        M_spikes_tmp = M_spikes;
                    end
                end 
                M_spikes_mean = mean(M_spikes, 2);     
                tmp = M_spikes_mean ./ bin_size;
                firing_rate = [firing_rate; tmp];
            end 
            pca_matrix = [pca_matrix, firing_rate]; 
        end 
        pca_matrix_array = [pca_matrix_array, pca_matrix]; 
    end
 
   
    % z-score
    X = zscore(pca_matrix_array, 0, 1); 
    
    % 10D FA
    is_finite_col = all(isfinite(X), 1);          % rimuovo colonne con NaN/Inf altrimenti factoran non funziona 
    std_col = std(X, 0, 1);                       % deviazione std per colonna
    th = 1e-8;                                    % soglia "costante"
    is_const_col = std_col <= th;                 % elimina neuroni a varianza ~0
    good_cols = is_finite_col & ~is_const_col;
    fprintf('FA: %d/%d canali esclusi (non finiti o std~0)\n', sum(~good_cols), numel(good_cols));

    X_full_cond{d}    = X;
    good_cols_cond{d} = good_cols; 

end

% Intersezione dei canali buoni nelle due condizioni 
good_cols_both = good_cols_cond{1} & good_cols_cond{2};
fprintf('Intersezione canali: %d totali (su %d)\n', sum(good_cols_both), size(X_full_cond{1},2));

% Ricalcolo X filtrate con la stessa maschera
X_motor = X_full_cond{1}(:, good_cols_both);
X_free  = X_full_cond{2}(:, good_cols_both);

% FA separata (10 fattori) 
k = 10;
[Lambda_m, Psi_m] = factoran(X_motor, k, 'Rotate','none', 'Maxit',1000);
[Lambda_f, Psi_f] = factoran(X_free,  k, 'Rotate','none', 'Maxit',1000);

% Concatenazione e SVD (Eq. 1) 
Wcat = [Lambda_m, Lambda_f];          % [N x 20] con lo stesso N in entrambe
[U,~,~] = svd(Wcat, 'econ');
U = U(:,1:2*k);                       % base ortonormale N x 20

    % Salva maschera e X "full" per ciascuna condizione
    if strcmp(ds_name, 'motor_BCI02.mat')
        good_cols_motor = good_cols;    
        X_motor_full    = X;            
    else
        good_cols_free  = good_cols;     
        X_free_full     = X;             
    end

    X_fa = X(:, good_cols);

    m = 10; 
    if strcmp(ds_name, 'motor_BCI02.mat')
        X_motor = X_fa; 
        [Lambda_m, Psi_m, ~, stats_m, F_m] = factoran(X_motor, 10, ...
            'Rotate','none', 'Scores','Bartlett', 'Maxit',1000);   
    else
        X_free = X_fa; 
        [Lambda_f, Psi_f, ~, stats_f, F_f] = factoran(X_free, 10, ...
            'Rotate','none', 'Scores','Bartlett', 'Maxit',1000);
    end
end 

%% Vediamo: (1) quanto bene lo spazio di una condizione spiega i propri dati; (2) quanto bene lo spazio di una condizione spiega i dati dell’altra.
% Ogni FA è stata stimata dopo aver scartato neuroni diversi; per confrontare spazi, bisogna usare lo stesso set di neuroni.
idx_full_common = find(good_cols_motor & good_cols_free);  % considero i neuroni buoni in entrambe le condizioni
Xm = X_motor_full(:, idx_full_common);                     % matrice motor che tiene in conto dei soli neuroni buoni in entrambe le condizioni
Xf = X_free_full(:,  idx_full_common);                     % matrice free che tiene in conto dei soli neuroni buoni in entrambe le condizioni

% Lambda_m e Lambda_f hanno righe solo per i neuroni buoni della rispettiva condizione.
rows_m_full = find(good_cols_motor);   % posizioni FULL dei neuroni usati in Lambda_m
rows_f_full = find(good_cols_free);    % posizioni FULL dei neuroni usati in Lambda_f

% pos_m e pos_f sono gli indici di riga per Lambda_m/Lambda_f che corrispondono esattamente alle colonne comuni di Xm/Xf.
[~, pos_m] = intersect(rows_m_full, idx_full_common, 'stable');
[~, pos_f] = intersect(rows_f_full, idx_full_common, 'stable');

% Sottoseleziona le loading delle due FA per gli stessi neuroni (n_common) in modo che le loading siano confrontabili.
Lam_m_c = Lambda_m(pos_m, :);   % (n_common x 10)
Lam_f_c = Lambda_f(pos_f, :);   % (n_common x 10)

% Ortonormalizza gli span e calcola %var spiegata
% Costruzione di due spazi fattoriali 
[Qm, ~] = qr(Lam_m_c, 0);
[Qf, ~] = qr(Lam_f_c, 0);

% Percentuale di varianza spiegata nello spazio 10D
% sum(var(Y,0,1): calcola la varianza di ogni colonna di Y e la somma per calcolare la varianza totale dei dati originali
% Y*Q: proietta i dati Y sui nuovi assi ortonormali (cambio di sistema di riferimento)
% sum(var(Y*Q,0,1)): calcolo della varianza dei dati proiettati sulla nuova base
pctvar = @(Y,Q) 100 * ( sum(var(Y*Q,0,1)) / sum(var(Y,0,1)) );

var_mm = pctvar(Xm, Qm); % quanta varianza di Motor è spiegata dai fattori Motor
var_mf = pctvar(Xf, Qm); % quanta varianza di Free-gaze è spiegata dai fattori Motor
var_ff = pctvar(Xf, Qf); % quanta varianza di Free-gaze è spiegata dai fattori Free-gaze
var_fm = pctvar(Xm, Qf); % quanta varianza di Motor è spiegata dai fattori Free-gaze

fprintf('\nVariance explained (%%):\n');
fprintf('Motor space on motor data:     %.1f%%\n', var_mm);
fprintf('Motor space on free-gaze data: %.1f%%\n', var_mf);
fprintf('Free space on free-gaze data:  %.1f%%\n', var_ff);
fprintf('Free space on motor data:      %.1f%%\n\n', var_fm);

%% Figure (1)
vals = [var_mm, var_fm;   % Motor data / Motor vs Free space
        var_mf, var_ff];  % Free data  / Motor vs Free space

figure('Color','w'); 
bar(vals);
ylabel('% variance explained');
clear set
set(gca, 'XTickLabel', {'Motor data','Free-gaze data'});
legend({'Motor FA space','Free-gaze FA space'});
title('Cross-condition FA variance explained');

%% Figure (2)
M = [var_mm, var_mf;   % fattori dalle righe, dati dalle colonne
     var_fm, var_ff];  %    (Motor, Free) x (Motor, Free)

figure('Color','w'); 
imagesc(M); axis equal tight
cmap = parula(256);
cmap = cmap(1:230, :); % rimuove i colori finali (più gialli)
colormap(cmap); c = colorbar; c.Label.String = '% variance explained';
xticks([1 2]); xticklabels({'Motor data','Free-gaze data'});
yticks([1 2]); yticklabels({'Motor FA space','Free-gaze FA space'});
title('Cross-condition FA variance explained');

for i=1:2, for j=1:2
    text(j,i,sprintf('%.1f%%',M(i,j)),'HorizontalAlignment','center',...
        'FontWeight','bold','Color','w');
end, end

%% Joint space
% Costruzione dello spazio congiunto dai due set di loading FA
% Normalizziamo le due matrici di loading dividendo per radice della somma dei 
% quadrati di tutti gli elementi, in modo da dare pari peso alle due condizioni.
Lm_n = Lam_m_c ./ norm(Lam_m_c, 'fro');
Lf_n = Lam_f_c ./ norm(Lam_f_c, 'fro');
L_joint = [Lambda_m, Lambda_f];                  

% SVD → base ortonormale delle "original dimensions"
[U_joint, S_joint, ~] = svd(L_joint, 'econ');

% Percentuale di varianza spiegata nello spazio 20D
pctvar = @(Y,U_joint) 100 * ( sum(var(Y*U_joint, 0, 1)) / sum(var(Y, 0, 1)) );

comb_motor = pctvar(Xm, U_joint);            % var di MOTOR spiegata dal 20D congiunto
comb_free  = pctvar(Xf, U_joint);            % var di FREE  spiegata dal 20D congiunto

fprintf('Combined 20D space variance explained (%%):\n');
fprintf('  Motor data: %.1f%%\n', comb_motor);
fprintf('  Free-gaze data: %.1f%%\n', comb_free);
U20 = U_joint; 

%% Figure
% Proiezioni per condizione
Z_motor = Xm * U20;          % (samples x 20)
Z_free  = Xf * U20;          % (samples x 20)
r = size(U20,2);

% Varianza per dimensione
v_m = var(Z_motor, 0, 1);    % 1 x 20
v_f = var(Z_free,  0, 1);    % 1 x 20

% Normalizza a percentuale per condizione (somma = 100% ciascuna)
exp_m = 100 * v_m / sum(v_m);
exp_f = 100 * v_f / sum(v_f);

% Prepara matrice per barplot raggruppato (20 coppie di barre)
K = min(20, size(U20,2));
Y = [exp_m(1:K); exp_f(1:K)].';   % K x 2 (col1: motor, col2: free)

figure('Color','w'); clf
hb = bar(1:K, Y, 'grouped'); hold on
hb(1).FaceColor = [0.20 0.50 1.00];   % blu (action/motor)
hb(2).FaceColor = [1.00 0.40 0.10];   % arancione (imagery/free-gaze)
hb(1).BarWidth  = 0.9; hb(2).BarWidth = 0.9;

ylim([0 70]); xlim([0.5 K+0.5])
xticks(1:K); xticklabels(string(1:K))
ylabel('% Variance'); xlabel('Dimensions')
legend({'motor-only','free-gaze'}, 'Location','northeast')
title('Variance explained in shared 20D space')
set(gca,'Box','off','FontSize',12)

%% free-unique, motor-unique, shared
% Covarianze nelle 20D 
Sm = cov(Z_motor, 1);                
Sf = cov(Z_free,  1);

Sm = Sm / trace(Sm);                 % bilancia le scale tra condizioni
Sf = Sf / trace(Sf);

% Autovettori generalizzati: Sf v = λ Sm v
[Vec, D] = eig(Sf, Sm);
lambda = real(diag(D));

% Ordina per "forza di separazione" (|log λ|)
[~, ord] = sort(abs(log(lambda)), 'descend');
V = real(Vec(:, ord));

% Ortonormalizza le direzioni nello spazio 20D
[Q_all, ~] = qr(V, 0);                % r x r, colonne ortonormali
z = log(lambda(ord));                 % simmetrico intorno a 0 se bilanciato

% --- Banda shared (più larga) attorno a 0 su log λ
eps_shared = 0.65;                     % aumenta (es. 0.5–0.6) per più shared
is_shared = abs(z) <= eps_shared;
is_free   =  z >  eps_shared;
is_motor  =  z < -eps_shared;

Q_shared = Q_all(:, is_shared);
Q_free   = Q_all(:, is_free);
Q_motor  = Q_all(:, is_motor);

fprintf('Subspace dims -> motor-unique: %d, shared: %d, free-unique: %d\n', ...
        size(Q_motor,2), size(Q_shared,2), size(Q_free,2));

% % var spiegata da ciascun sottospazio (per condizione)
pctvar = @(Y,Q) 100 * ( sum(var(Y*Q,0,1)) / sum(var(Y,0,1)) );

mu_motor = pctvar(Z_motor, Q_motor);
sh_motor = pctvar(Z_motor, Q_shared);
fu_motor = pctvar(Z_motor, Q_free);

mu_free  = pctvar(Z_free,  Q_motor);
sh_free  = pctvar(Z_free,  Q_shared);
fu_free  = pctvar(Z_free,  Q_free);

fprintf('\n%% var explained by subspaces:\n');
fprintf('  Motor data   -> Mu: %.1f, Sh: %.1f, Fu: %.1f\n', mu_motor, sh_motor, fu_motor);
fprintf('  Free-gaze    -> Mu: %.1f, Sh: %.1f, Fu: %.1f\n\n', mu_free,  sh_free,  fu_free);

%% Figure (1)
% CONCATENA varianze in ordine: [motor-unique  free-unique  shared]
motor_vals = [mu_motor, fu_motor, sh_motor];
free_vals  = [mu_free,  fu_free,  sh_free];

% dimensioni per ciascun blocco
n_mu = length(mu_motor);
n_fu = length(fu_motor);
n_sh = length(sh_motor);

x = 1:length(motor_vals);   % asse x per le barre

figure('Color','w'); clf; hold on;

% barre affiancate stile paper
bar(x-0.15, motor_vals, 0.3, 'FaceColor', [0.2 0.4 1]);  % blu   (motor dataset)
bar(x+0.15, free_vals,  0.3, 'FaceColor', [1 0.4 0.1]);  % arancione (free dataset)

ylabel('% variance');
xlabel('Subspace dimensions');
title('% variance per dimension — motor-unique / free-unique / shared');

ylim([0 max([motor_vals free_vals])*1.2]);
xlim([0 length(x)+1]);
set(gca,'TickDir','out'); box off;

% separatori visivi dei blocchi
x_mu_end = n_mu;
x_fu_end = n_mu + n_fu;

line([x_mu_end+0.5 x_mu_end+0.5], ylim, 'Color',[0.7 0.7 0.7],'LineStyle','--');
line([x_fu_end+0.5 x_fu_end+0.5], ylim, 'Color',[0.7 0.7 0.7],'LineStyle','--');

% etichette categorie sotto
xticks([mean(1:x_mu_end), mean(x_mu_end+1:x_fu_end), mean(x_fu_end+1:length(x))]);
xticklabels({'motor-unique','free-unique','shared'});

legend({'Motor data','Free-gaze data'}, ...
       'Location','northoutside','Orientation','horizontal');
