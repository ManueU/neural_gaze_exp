% =========================================================
% DESCRIZIONE:
% Applica una Canonical Correlation Analysis (CCA) per valutare:
% - quanto le traiettorie neurali di due condizioni sperimentali 
%   condividano una struttura dinamica comune nello stesso spazio latente;
% - la similarità tra le condizioni tramite le correlazioni canoniche;
% - la forma delle traiettorie nei primi assi canonici (Can1–Can3).
%
% La procedura prevede:
% 1) Z-score delle condizioni rispetto ai parametri di normalizzazione (muZ, sigmaZ);
% 2) proiezione dei dati sulle componenti principali comuni (coeff);
% 3) estrazione delle prime PC (pc_idx_cca) come input per la CCA;
% 4) calcolo di U e V, le traiettorie canoniche per le due condizioni;
% 5) visualizzazione tridimensionale delle traiettorie smussate e 
%    rappresentazione delle correlazioni canoniche.
%
% Prima dell’esecuzione, assicurarsi di impostare correttamente:
%   condition_sep : cell array con le due condizioni da confrontare (run codice pca_w_fixed.m)
%   muZ, sigmaZ   : parametri di normalizzazione (run codice pca_w_fixed.m)
%   coeff         : componenti principali comuni (PCA condivisa) (run codice pca_w_fixed.m)
%   pc_idx_cca    : indici delle PC da includere nella CCA
%   n_targets     : numero di target comportamentali
%   nbin, w       : parametri per segmentazione e smoothing
% =========================================================


%% CCA alignment
Xz_cond1 = (condition_sep{1} - muZ) ./ sigmaZ;
Xz_cond2 = (condition_sep{2} - muZ) ./ sigmaZ;

% Proiezione sulle PC comuni
scores1 = Xz_cond1 * coeff;   
scores2 = Xz_cond2 * coeff;

pc_idx_cca = 1:10;

Xcca = scores1(:, pc_idx_cca);   % condizione 1
Ycca = scores2(:, pc_idx_cca);   % condizione 2
Xcca = zscore(Xcca);
Ycca = zscore(Ycca);

[A, B, r, U, V, stats] = canoncorr(Xcca, Ycca);
disp('Canonical correlations (CCA):');
disp(r);


%% Figure (1) - trajectories in canonical space)
figure('Color','w'); hold on;
axis equal; grid on;
xlabel('Can1'); ylabel('Can2'); zlabel('Can3');
title('Trajectories in canonical space (CCA)');

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
end

h1 = plot3(nan,nan,nan,'-k','LineWidth',2,'DisplayName','Condizione 1');
h2 = plot3(nan,nan,nan,'--k','LineWidth',2,'DisplayName','Condizione 2');
legend([h1 h2],'Location','northeastoutside');


%% Figure (2)
figure('Color','w');
bar(r, 'FaceColor', [0.6 0.8 0.6], 'EdgeColor','none');
ylim([0 1])
xlabel('Canonical dimension');
ylabel('Correlation');
title('Canonical correlations between conditions');
box("off")


%% Figure (3)
r_thr = 0.5; 
figure('Color','w');
bar(r, 'FaceColor', [0.6 0.8 0.6], 'EdgeColor','none');
yline(r_thr, '--', 'Shared threshold');
ylim([0 1])
xlabel('Canonical dimension');
ylabel('Correlation');
title('Canonical correlations between conditions');
box("off")
