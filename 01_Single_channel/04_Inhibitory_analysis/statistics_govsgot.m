%% Scelta automatica del test paired in base alla normalità delle differenze
% Richiede nel workspace:
%   common_by_target
%
% Per ogni target:
%   1) prende i soli canali comuni tra GO e GOT
%   2) calcola le differenze d = GO - GOT
%   3) testa la normalità di d
%   4) se normale -> paired t-test
%      altrimenti -> Wilcoxon signed-rank test

clearvars -except common_by_target

measureVar = 'peak_inh_amp';
alphaNorm  = 0.05;   % soglia per test di normalità
alphaTest  = 0.05;   % soglia per test statistico

fn = fieldnames(common_by_target);
nTargets = numel(fn);

% Preallocazione
target_id      = cell(nTargets,1);
n_channels     = nan(nTargets,1);

mean_go        = nan(nTargets,1);
mean_got       = nan(nTargets,1);
mean_diff      = nan(nTargets,1);
median_diff    = nan(nTargets,1);

normality_h    = nan(nTargets,1);   % 0 = non rifiuto normalità, 1 = rifiuto
normality_p    = nan(nTargets,1);

selected_test  = strings(nTargets,1);
test_h         = nan(nTargets,1);
test_p         = nan(nTargets,1);
test_stat      = nan(nTargets,1);   % t-stat oppure signed-rank statistic
df_or_z        = nan(nTargets,1);   % df per t-test, z per signrank se disponibile

for i = 1:nTargets
    S = common_by_target.(fn{i});

    Tgo  = S.go;
    Tgot = S.got;

    target_id{i} = char(string(S.target_id));

    % Controllo presenza variabile
    if ~ismember(measureVar, Tgo.Properties.VariableNames) || ...
       ~ismember(measureVar, Tgot.Properties.VariableNames)
        warning('Target %s: variabile %s non trovata. Salto.', ...
            string(S.target_id), measureVar);
        continue
    end

    % Controllo stesso numero di righe
    if height(Tgo) ~= height(Tgot)
        warning('Target %s: numero di righe diverso tra GO e GOT. Salto.', ...
            string(S.target_id));
        continue
    end

    % Controllo allineamento canali
    varsToCheck = intersect({'array_name','channel_global','channel_local'}, ...
                            Tgo.Properties.VariableNames);

    aligned = true;
    for v = 1:numel(varsToCheck)
        vn = varsToCheck{v};
        if ~isequal(Tgo.(vn), Tgot.(vn))
            warning('Target %s: GO e GOT non allineate per %s. Salto.', ...
                string(S.target_id), vn);
            aligned = false;
            break
        end
    end
    if ~aligned
        continue
    end

    % Estrai dati paired
    x = Tgo.(measureVar);
    y = Tgot.(measureVar);

    % Tieni solo coppie valide
    valid = ~isnan(x) & ~isnan(y);
    x = x(valid);
    y = y(valid);

    n = numel(x);
    n_channels(i) = n;

    if n < 3
        warning('Target %s: meno di 3 coppie valide. Test non eseguito.', ...
            string(S.target_id));
        continue
    end

    % Statistiche descrittive
    d = x - y;
    mean_go(i)     = mean(x);
    mean_got(i)    = mean(y);
    mean_diff(i)   = mean(d);
    median_diff(i) = median(d);

    %% Test di normalità sulle differenze
    % H = 0 -> normalità non rifiutata
    % H = 1 -> normalità rifiutata
    try
        [h_norm, p_norm] = lillietest(d, 'Alpha', alphaNorm);
    catch
        % fallback se lillietest non fosse disponibile
        [h_norm, p_norm] = jbtest(d, alphaNorm);
    end

    normality_h(i) = h_norm;
    normality_p(i) = p_norm;

    %% Scelta del test
    if h_norm == 0
        % Differenze compatibili con normalità -> paired t-test
        [h_t, p_t, ~, stats_t] = ttest(x, y, 'Alpha', alphaTest);

        selected_test(i) = "paired_ttest";
        test_h(i)        = h_t;
        test_p(i)        = p_t;
        test_stat(i)     = stats_t.tstat;
        df_or_z(i)       = stats_t.df;

    else
        % Differenze non normali -> Wilcoxon signed-rank
        [p_sr, h_sr, stats_sr] = signrank(x, y, 'Alpha', alphaTest);

        selected_test(i) = "signrank";
        test_h(i)        = h_sr;
        test_p(i)        = p_sr;

        if isfield(stats_sr, 'signedrank')
            test_stat(i) = stats_sr.signedrank;
        end
        if isfield(stats_sr, 'zval')
            df_or_z(i) = stats_sr.zval;
        end
    end
end

%% Tabella risultati
results_paired = table( ...
    target_id, ...
    n_channels, ...
    mean_go, ...
    mean_got, ...
    mean_diff, ...
    median_diff, ...
    normality_h, ...
    normality_p, ...
    selected_test, ...
    test_h, ...
    test_p, ...
    test_stat, ...
    df_or_z);

disp(results_paired)