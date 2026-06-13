%% Mann–Whitney (ranksum) GO vs GOT per ciascun target

clearvars -except T_go T_got

measureVar = 'peak_inh_amp';
alpha = 0.05;

%% --- Setup ---
T_go.array_name  = string(T_go.array_name);
T_got.array_name = string(T_got.array_name);

if iscell(T_go.target_id),  T_go.target_id  = string(T_go.target_id);  end
if iscell(T_got.target_id), T_got.target_id = string(T_got.target_id); end

targets = union(unique(T_go.target_id), unique(T_got.target_id), 'stable');
nTargets = numel(targets);

%% --- Preallocazione ---
target_id   = cell(nTargets,1);
n_GO        = nan(nTargets,1);
n_GOT       = nan(nTargets,1);

median_GO   = nan(nTargets,1);
median_GOT  = nan(nTargets,1);

p_value     = nan(nTargets,1);
h_value     = nan(nTargets,1);
z_value     = nan(nTargets,1);

effect_rbc  = nan(nTargets,1); % rank-biserial correlation

%% --- Loop per target ---
for i = 1:nTargets
    tgt = targets(i);
    target_id{i} = char(tgt);

    % Estrai dati
    goVals  = T_go{T_go.target_id == tgt, measureVar};
    gotVals = T_got{T_got.target_id == tgt, measureVar};

    % Rimuovi NaN
    goVals  = goVals(~isnan(goVals));
    gotVals = gotVals(~isnan(gotVals));

    n_GO(i)  = numel(goVals);
    n_GOT(i) = numel(gotVals);

    if n_GO(i) < 2 || n_GOT(i) < 2
        warning('Target %s: dati insufficienti', tgt)
        continue
    end

    % Mediane
    median_GO(i)  = median(goVals);
    median_GOT(i) = median(gotVals);

    %% Mann–Whitney
    [p, h, stats] = ranksum(goVals, gotVals, 'alpha', alpha);

    p_value(i) = p;
    h_value(i) = h;

    if isfield(stats,'zval')
        z_value(i) = stats.zval;
    end

    %% Effect size: rank-biserial correlation
    % r = z / sqrt(N)
    N = n_GO(i) + n_GOT(i);
    effect_rbc(i) = z_value(i) / sqrt(N);
end

%% --- Correzione multipli (Bonferroni) ---
p_bonf = min(p_value * nTargets, 1);

%% --- Tabella risultati ---
results_mw = table( ...
    target_id, ...
    n_GO, n_GOT, ...
    median_GO, median_GOT, ...
    p_value, p_bonf, ...
    h_value, z_value, effect_rbc);

disp(results_mw)