% =========================================================
% DESCRIZIONE:
% Questo script:
% 1) Valuta la stabilità delle correlazioni canoniche (CCA) al variare
%    del numero di componenti principali m usate come input;
% 2) Esegue una CCA fissando m = 10 e ottiene le traiettorie canoniche
%    U e V per le due condizioni;
% 3) Visualizza le traiettorie canoniche (Can1–Can3) per ciascun target;
% 4) Mostra le correlazioni canoniche come grafico a barre, con e senza
%    soglia di "condivisione" (r_thr).
%
% La procedura prevede:
% - normalizzazione delle due condizioni (condition_sep{1}, condition_sep{2})
%   usando i parametri muZ e sigmaZ;
% - proiezione dei dati neurali sulle componenti principali comuni (coeff);
% - per ciascun valore di m in m_values:
%       * selezione delle prime m PC;
%       * bootstrap sulle righe (ricampionamento con rimpiazzamento);
%       * calcolo delle Canonical Correlations r;
%       * stima di media e deviazione standard delle prime nCC;
% - scelta di un valore di m (qui m = 10) e calcolo della CCA finale;
% - utilizzo di U e V per tracciare, per ogni target, le traiettorie
%   nelle prime tre dimensioni canoniche (Can1–Can3), con smoothing;
% - costruzione di grafici riassuntivi:
%       * errore standard delle prime nCC CC in funzione di m;
%       * barre delle correlazioni canoniche;
%       * barre con soglia r_thr per identificare dimensioni “condivise”.
%
% Prima dell’esecuzione, assicurarsi di impostare correttamente:
%   condition_sep : cell array con le due condizioni da confrontare
%                   (run codice pca_w_fixed.m)
%   muZ, sigmaZ   : parametri di normalizzazione (run pca_w_fixed.m)
%   coeff         : componenti principali comuni (PCA condivisa)
%                   (run pca_w_fixed.m)
%   m_values      : lista dei numeri di PC da testare nella CCA (stability)
%   nCC           : numero di canonical correlations da considerare (es. 4)
%   nBoot         : numero di bootstrap per la stima di stabilità
%   n_targets     : numero di target
%   nbin          : numero di bin temporali per target
%   w             : finestra per lo smoothing (smoothdata)
%   r_thr         : soglia sulle correlazioni canoniche per Figure (4)
% =========================================================


%% CCA stability
Xz_cond1 = (condition_sep{1} - muZ) ./ sigmaZ;
Xz_cond2 = (condition_sep{2} - muZ) ./ sigmaZ;

% Proiezione sulle PC comuni
scores1 = Xz_cond1 * coeff;   
scores2 = Xz_cond2 * coeff;

m_values = [4 6 8 10 12 15 20];
nCC = 4;
nBoot = 20;
r_all = nan(numel(m_values), nCC, nBoot);

for im = 1:numel(m_values)
    m = m_values(im);
    pc_idx_cca = 1:m;

    % Dati proiettati sulle prime m PC
    Xcca_full = scores1(:, pc_idx_cca);
    Ycca_full = scores2(:, pc_idx_cca);

    % Bootstrap sulla dimensione dei campioni (righe)
    nSamples = size(Xcca_full, 1);

    for b = 1:nBoot
        % Campionamento con rimpiazzamento delle righe
        idx = randsample(nSamples, nSamples, true);

        Xcca = zscore(Xcca_full(idx, :));
        Ycca = zscore(Ycca_full(idx, :));

        [A, B, r, U, V, stats] = canoncorr(Xcca, Ycca);

        % Salvo solo le prime nCC canonical correlations
        r_all(im, :, b) = r(1:nCC);
    end
end

% Calcolo media e deviazione standard delle prime 4 CC per ogni m
r_mean = mean(r_all, 3, 'omitnan');
r_std  = std(r_all,  0, 3, 'omitnan');

fprintf('Canonical correlations (prime %d CC) per diversi m:\n', nCC);
for im = 1:numel(m_values)
    m = m_values(im);
    fprintf('m = %2d:  ', m);
    for k = 1:nCC
        fprintf('CC%d = %.3f ± %.3f   ', k, r_mean(im,k), r_std(im,k));
    end
    fprintf('\n');
end


%% Figure (1)
figure('Color','White'); 
hold on;
for k = 1:nCC
    errorbar(m_values, r_mean(:,k), r_std(:,k), '-o', 'DisplayName', sprintf('CC%d',k));
end
ylim([0 1])
xlabel('PCs');
ylabel('Canonical correlation');
title(sprintf('Stability of %d canonical correlations', nCC));
legend('Location','best');
grid on;

%% CCA alignment with m = 10
pc_idx_cca = 1:10;

Xcca = scores1(:, pc_idx_cca);   % condizione 1
Ycca = scores2(:, pc_idx_cca);   % condizione 2
Xcca = zscore(Xcca);
Ycca = zscore(Ycca);

[A, B, r, U, V, stats] = canoncorr(Xcca, Ycca);
disp('Canonical correlations (CCA):');
disp(r);

% --------- CCA-distance per la CCA finale (m = 10) ----------------------
% Usiamo le prime k_F CC (qui k_F = min(nCC, length(r)))
k_F = min(nCC, numel(r));
d_CCA_final = 1 - mean(r(1:k_F));

fprintf('\nCCA-distance (prime %d CC, m = 10): %.4f\n', k_F, d_CCA_final);
% ------------------------------------------------------------------------

%% Figure (2) - trajectories in canonical space
figure('Color','w'); 
hold on;
axis equal; 
grid on;
xlabel('Can1'); ylabel('Can2'); zlabel('Can3');
title('CCA');

for tg = 1:n_targets
    idx = (tg-1)*nbin + (1:nbin);

    Ut = U(idx, 1:3);   % cond1 in spazio canonico
    Vt = V(idx, 1:3);   % cond2 in spazio canonico

    % smoothing
    u1 = smoothdata(Ut(:,1), 'gaussian', w);
    u2 = smoothdata(Ut(:,2), 'gaussian', w);
    u3 = smoothdata(Ut(:,3), 'gaussian', w);

    v1 = smoothdata(Vt(:,1), 'gaussian', w);
    v2 = smoothdata(Vt(:,2), 'gaussian', w);
    v3 = smoothdata(Vt(:,3), 'gaussian', w);

    % condizione 1 (linea continua)
    plot3(u1, u2, u3, '-', ...
        'Color', colors(tg,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('%d - cond1', tg));

    % condizione 2 (linea tratteggiata)
    plot3(v1, v2, v3, '--', ...
        'Color', colors(tg,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('%d - cond2', tg));

    % marker inizio/fine per condizione 1
    plot3(u1(1), u2(1), u3(1), 'o', ...
              'MarkerSize', 8, 'MarkerFaceColor', colors(tg,:), ...
              'MarkerEdgeColor', colors(tg,:), 'HandleVisibility','off');
    plot3(u1(end), u2(end), u3(end), '^', ...
              'MarkerSize', 8, 'MarkerFaceColor', colors(tg,:), ...
              'MarkerEdgeColor', colors(tg,:), 'HandleVisibility','off');

    % marker inizio/fine per condizione 2 allineata
    plot3(v1(1), v2(1), v3(1), 'o', ...
              'MarkerSize', 8, 'MarkerFaceColor', colors(tg,:), ...
              'MarkerEdgeColor', colors(tg,:), 'HandleVisibility','off');
    plot3(v1(end), v2(end), v3(end), '^', ...
              'MarkerSize', 8, 'MarkerFaceColor', colors(tg,:), ...
              'MarkerEdgeColor', colors(tg,:), 'HandleVisibility','off');
end

hStart = plot3(nan,nan,nan, 'o', 'Color','k', 'MarkerFaceColor','k');
hEnd   = plot3(nan,nan,nan, '^', 'Color','k', 'MarkerFaceColor','k');

styles = {'-','--'};
condNames = cell(1,2);
for c = 1:2
    [~, baseName] = fileparts(filename{c});
    [namePart, ~] = strtok(baseName, '_');
    condNames{c} = namePart;
    ls = styles{c};
    hStyles(c) = plot3(nan,nan,nan, ls, 'Color','k', 'LineWidth', 2);
end
% 
% legend([hStyles hStart hEnd], {condNames{:}, 'Start', 'End'}, ...
%        'Location','northeastoutside');


%% Figure (3)
figure('Color','w');
bar(r, 'FaceColor', [0.6 0.8 0.6], 'EdgeColor','none');
ylim([0 1])
xlabel('Canonical dimension');
ylabel('Correlation');
title('Canonical correlations between conditions');
box("off")


%% Figure (4)
r_thr = 0.5; 
figure('Color','w');
bar(r, 'FaceColor', [0.6 0.8 0.6], 'EdgeColor','none');
yline(r_thr, '--', 'Shared threshold');
ylim([0 1])
xlabel('Canonical dimension');
ylabel('Correlation');
title('Canonical correlations between conditions');
box("off")


