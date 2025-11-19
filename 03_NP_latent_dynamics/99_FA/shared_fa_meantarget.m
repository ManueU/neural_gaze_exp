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
    X_full_raw{d} =  pca_matrix_array; 
    X = zscore(pca_matrix_array, 0, 1); 
    X_full_cond{d}    = X;


    is_finite_col = all(isfinite(X), 1);          % rimuovo colonne con NaN/Inf altrimenti factoran non funziona 
    std_col = std(X, 0, 1);                       % deviazione std per colonna
    th = 1e-8;                                    % soglia "costante"
    is_const_col = std_col <= th;                 % elimina neuroni a varianza ~0
    good_cols = is_finite_col & ~is_const_col;
    good_cols_cond{d} = good_cols; 
    fprintf('FA: %d/%d canali esclusi (non finiti o std~0)\n', sum(~good_cols), numel(good_cols));

end

% Intersezione dei canali buoni nelle due condizioni 
good_cols_both = good_cols_cond{1} & good_cols_cond{2};
fprintf('Intersezione canali: %d totali (su %d)\n', sum(good_cols_both), size(X_full_cond{1},2));

% Ricalcolo X filtrate con la stessa maschera
X_motor_raw = X_full_raw{1}(:, good_cols_both);
X_free_raw = X_full_raw{2}(:, good_cols_both);
X_motor = X_full_cond{1}(:, good_cols_both);
X_free  = X_full_cond{2}(:, good_cols_both);

k = 10;
[Lambda_m, Psi_m] = factoran(X_motor, k, 'Rotate','none', 'Maxit',1000);
[Lambda_f, Psi_f] = factoran(X_free,  k, 'Rotate','none', 'Maxit',1000);

%% Joint space
% Concatenazione e SVD (Eq. 1) 
Wcat = [Lambda_m, Lambda_f];            % [N x 20] con lo stesso N in entrambe
[U,~,~] = svd(Wcat, 'econ');
U = U(:,1:2*k);  
[Qm,~] = qr(Lambda_m,0);                % N x 10 ortonormale
[Qf,~] = qr(Lambda_f,0);                % N x 10 ortonormale

W_motor = Qm;                           % N x 10
W_free  = Qf;                           % N x 10
W_comb  = U;                            % N x 20  (dalla SVD combinata)

%% Vediamo: 
% (1) quanto bene lo spazio di una condizione spiega i propri dati; 
% (2) quanto bene lo spazio di una condizione spiega i dati dell’altra.

% mean-centering
Z_motor = X_motor_raw - mean(X_motor_raw,1); 
Z_free = X_free_raw - mean(X_free_raw,1); 

varExpl = @(Z, W) trace((Z*W)'*(Z*W)) / trace(Z'*Z) * 100;
p_motor_inMotor = varExpl(Z_motor, W_motor);   % var. motor in proprio spazio
p_free_inFree   = varExpl(Z_free,  W_free);    % var. free in proprio spazio

p_motor_inFree  = varExpl(Z_motor, W_free);    % cross: motor nello spazio free
p_free_inMotor  = varExpl(Z_free,  W_motor);   % cross: free nello spazio motor

p_motor_inComb = varExpl(Z_motor, U);
p_free_inComb  = varExpl(Z_free,  U);

fprintf('\nVarianza spiegata (10D ans 20D FA)\n');
fprintf('Free 10D : free %.1f%% | motor %.1f%%\n', p_free_inFree,  p_motor_inFree);
fprintf('Motor   10D : free %.1f%% | motor %.1f%%\n', p_free_inMotor, p_motor_inMotor);
fprintf('Comb 20D: free %.1f%% | motor %.1f%%\n', p_free_inComb,  p_motor_inComb);


%% Figure (1)
vals = [p_motor_inMotor, p_motor_inFree;   % Motor data / Motor vs Free space
        p_free_inMotor, p_free_inFree];  % Free data  / Motor vs Free space

figure('Color','w'); 
bar(vals);
ylabel('% variance explained');
clear set
set(gca, 'XTickLabel', {'Motor data','Free-gaze data'});
legend({'Motor FA space','Free-gaze FA space'});
title('Cross-condition FA variance explained');

%% Figure (2)
M = [p_motor_inMotor, p_free_inMotor;   % fattori dalle righe, dati dalle colonne
     p_motor_inFree, p_free_inFree];    %    (Motor, Free) x (Motor, Free)

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

%% Prijection of data into the 20D space 
% (1) normalizzazione per canale come nel paper (somma var tra condizioni)
v_motor = var(X_motor_raw,0,1);
v_free  = var(X_free_raw, 0,1);
v_sum   = v_motor + v_free;

Z_motor_n = X_motor_raw ./ sqrt(v_sum);
Z_free_n  = X_free_raw  ./ sqrt(v_sum);

% (2) proiezione nello spazio combinato (20D)
Z_motor_proj = Z_motor_n * U;   % T x 20
Z_free_proj  = Z_free_n  * U;   % T x 20

% (3) % di varianza spiegata PER DIMENSIONE (Eq.2, senza sommare)
totVar_motor = trace(Z_motor_n' * Z_motor_n);
totVar_free  = trace(Z_free_n'  * Z_free_n);

perc_motor = diag(Z_motor_proj' * Z_motor_proj) ./ totVar_motor * 100;    % 20x1
perc_free  = diag(Z_free_proj'  * Z_free_proj)  ./ totVar_free  * 100;   % 20x1


%% Figure (3)
figure('Color','w'); 
b = bar([perc_motor(:), perc_free(:)], 'grouped'); % 20 gruppi (dimensioni)
b(1).FaceColor = [0.12 0.47 0.71];  % blu
b(2).FaceColor = [0.89 0.39 0.13];  % arancione
b(1).EdgeColor = 'none'; b(2).EdgeColor = 'none';

% assi/etichette
xlabel('original dimensions');
ylabel('% variance');
legend({'motor activity','free activity'}, 'Location','northeast');
clear set
set(gca,'Box','off','TickDir','out','FontSize',11);

xlim([0.5 20.5]);
yl = max([perc_motor; perc_free]);
ylim([0, max(40, ceil(yl/5)*5)]);   % auto; minimo 40% come nel pannello

for ii=1:numel(b)
    b(ii).BarWidth = 0.9;
end


%% free-unique, motor-unique, shared
% Covarianze nelle 20D 
[T1,d] = size(Z_motor_proj); [T2,d2] = size(Z_free_proj);  assert(d==d2);

C1 = cov(Z_motor_proj,1);  
C2 = cov(Z_free_proj,1);  

alpha_motor = 0.03;                 % PIÙ permissiva quando cerchi v1
alpha_free  = 0.0075;               % PIÙ stretta quando cerchi v2

nu1 = alpha_motor * trace(C2);      % vincolo su C2 per v1 (motor-unique)
nu2 = alpha_free  * trace(C1);      % vincolo su C1 per v2 (free-unique)

X1c = Z_motor_proj;  X2c = Z_free_proj;  dcur = d;
Q1 = []; Q2 = [];       
var1 = []; var2 = [];

max_comp = d;     
found = 0;
while true
    disp(found)
    % Covarianze nello spazio residuo
    C1c = cov(X1c,1);  
    C2c = cov(X2c,1);  

    % Obiettivo: trovare un vettore v che 1) massimizzi la varianza del
    % task1 (v'C1v grande), 2) non abbiamo troppa varianza nel task 2 (v'C2v <= soglia).

    [eCo_min_vec, eCo_min_val] = eigs(C2c, 1, 'smallestreal');              % Ricerco gli autovalori reali più piccoli di C2c
    v1 = []; v2 = [];
    if eCo_min_val <= nu1                                                    % Verifica necessaria perché: se anche la varianza minima possibile di C2 è più grande della soglia,..
        mu_lo = 0; mu_hi = 1e6; v1_ok = []; valCo1_ok = Inf;                % allora nessuna direzione v può rispettare il vincolo (2)

        [v_tmp,~] = eigs(C1c, 1, 'largestreal');                            % Ricerco l'autovettore principale di C1c, quello associato al più grande autovalore reale
        v_tmp = v_tmp/norm(v_tmp);
        if real(v_tmp' * C2c * v_tmp) <= nu1                                 % Verifico nuovamente il vincolo su C2c (2): quanta varianza task2 ha nella stessa direzione?
            v1 = v_tmp;                                                     % Se la varianza è già sotto la soglia, ho trovato il vettore che soddisfi (1) e (2) (vettore unique per il task1)
        else
            % Risolvo il problema di ottimizzazione utilizzando il metodo
            % dei moltiplicatori di Lagrange: bisogna trovare un
            % compromesso tra massimizzare la varianza di task1 e
            % minimizare quella di task2. Il problema si sposta sulla
            % ricerca del mu che bislanci i due obiettivi. 
            for it=1:60
                mu = 0.5*(mu_lo+mu_hi);
                [v_mu,~] = eigs(C1c - mu*C2c, 1, 'largestreal');             
                v_mu = v_mu / norm(v_mu);
                valCo = real(v_mu' * C2c * v_mu);
                if valCo <= nu1
                    v1_ok = v_mu; valCo1_ok = valCo; mu_hi = mu;            % Calcolo la varianza nel task2 in direzione trovata e mi salvo i vettori che soddisfano la condizione (2)...
                else                                                        % Riduco mu in modo da trovare una soluzione con ancora più varianza in C1c.
                    mu_lo = mu;
                end
                if mu_hi - mu_lo < 1e-8, break; end
            end
            if ~isempty(v1_ok), v1 = v1_ok; end                             % Salvo il vettore v1 che spiega tanta varianza nel task1 e poca varianza nel task2.
        end
    end

    % Ripeto la stessa operazione per trovare il vettore v2 che spiega
    % tanta varianza nel task2 e poca varianza nel task1.
    [eC1_min_vec, eC1_min_val] = eigs(C1c, 1, 'smallestreal');
    if eC1_min_val <= nu2
        mu_lo = 0; mu_hi = 1e6; v2_ok = []; valCo2_ok = Inf;

        [v_tmp,~] = eigs(C2c, 1, 'largestreal'); 
        v_tmp = v_tmp/norm(v_tmp);
        if real(v_tmp' * C1c * v_tmp) <= nu2
            v2 = v_tmp;
        else
            for it=1:60
                mu = 0.5*(mu_lo+mu_hi);
                [v_mu,~] = eigs(C2c - mu*C1c, 1, 'largestreal');
                v_mu = v_mu / norm(v_mu);
                valCo = real(v_mu' * C1c * v_mu);
                if valCo <= nu2
                    v2_ok = v_mu; valCo2_ok = valCo; mu_hi = mu;
                else
                    mu_lo = mu;
                end
                if mu_hi - mu_lo < 1e-8, break; end
            end
            if ~isempty(v2_ok), v2 = v2_ok; end
        end
    end

    % Verifico se i vettori trovati sono validi 
    ok1 = ~isempty(v1) && (real(v1' * C2c * v1) <= nu1);
    ok2 = ~isempty(v2) && (real(v2' * C1c * v2) <= nu2);
    if ~ok1 && ~ok2
        break;  % Stop: niente altri vettori unique che rispettino il vincolo
    end

    % Se entrambi i vettori rispettano i vincoli, scegliamo quello che spiega più varianza nel proprio task
    if ok1 && ok2
        pick1 = real(v1' * C1c * v1);
        pick2 = real(v2' * C2c * v2);
        pick_task = (pick1 >= pick2) + 1;   % 1->task1, 2->task2
    elseif ok1
        pick_task = 1;
    else
        pick_task = 2;
    end

    if pick_task == 1
        Q1 = [Q1, v1]; 
        var1 = [var1; real(v1' * C1c * v1)];
        n  = length(v1); 
        P = eye(n) - v1*v1';    % deflazione ortogonale
    else
        Q2 = [Q2, v2]; 
        var2 = [var2; real(v2' * C2c * v2)];
        n  = length(v2); 
        P = eye(n) - v2*v2';    % deflazione ortogonale
    end

    % proietta i dati nel sottospazio ortogonale (togli la dimensione scelta)
    X1c = X1c * P;
    X2c = X2c * P;

    found = found + 1;
    if found >= max_comp
        break;                  % fermati comunque dopo d deflazioni
    end
end

% Base dello spazio shared (complemento ortogonale di [Q1 Q2])
Qall = [Q1, Q2];
if isempty(Qall)
    Q_shared = eye(d);
else
    Psh = eye(size(Qall,1)) - Qall*(Qall');   % proiettore sul complementare
    [Q_shared,~] = qr(Psh, 0);                % base ortonormale dello shared
    if isempty(Q_shared)                      % protezione degenerate
        Q_shared = orth(Psh);
    end
end


%% Varianza spiegata motor-unique, free-unique, shared
d = size(Z_motor_proj,2);
I = eye(d);

Q1o = orth(Q1);                    % 1) ortonormalizza Q1
P1  = I - Q1o*Q1o.';

Q2o = orth(P1 * Q2);               % 2) ortogonalizza Q2 rispetto a Q1
P12 = I - Q1o*Q1o.' - Q2o*Q2o.';

Qsho = null([Q1o Q2o].');       % 3) shared ortogonale a (Q1,Q2)

P1  = Q1o * Q1o.';      % proiettori veri (idempotenti e simmetrici)
P2  = Q2o * Q2o.';
Psh = Qsho* Qsho.';

den_m = trace(Z_motor_proj' * Z_motor_proj);
den_f = trace(Z_free_proj'  * Z_free_proj);

motor_u1 = 100 * trace(Z_motor_proj * P1  * Z_motor_proj') / den_m;
motor_u2 = 100 * trace(Z_motor_proj * P2  * Z_motor_proj') / den_m;
motor_sh = 100 * trace(Z_motor_proj * Psh * Z_motor_proj') / den_m;

free_u1  = 100 * trace(Z_free_proj  * P1  * Z_free_proj' ) / den_f;
free_u2  = 100 * trace(Z_free_proj  * P2  * Z_free_proj' ) / den_f;
free_sh  = 100 * trace(Z_free_proj  * Psh * Z_free_proj' ) / den_f;

fprintf('Somma Motor = %.1f\n', motor_u1+motor_u2+motor_sh);
fprintf('Somma Free  = %.1f\n', free_u1 +free_u2 +free_sh);


%% Figure (4)
figure('Color', 'w');
vals_m = [motor_u1, motor_u2, motor_sh];
vals_f = [free_u1,  free_u2,  free_sh];

barData = [vals_m; vals_f];
b = bar(barData, 'stacked');
b(1).FaceColor = [0.2 0.6 1.0];   % motor-unique
b(2).FaceColor = [1.0 0.4 0.4];   % free-unique
b(3).FaceColor = [0.6 0.6 0.6];   % shared

clear set
set(gca,'XTickLabel',{'Motor task','Free-gaze task'}, 'FontSize',12);
ylabel('% Variance explained');
legend({'Motor-unique','Free-unique','Shared'}, 'Location','northoutside','Orientation','horizontal');
ylim([0 100]);
title('Variance explained by each subspace');



