% =========================================================
% DESCRIZIONE:
% Applica un allineamento Procrustes per confrontare la forma delle 
% traiettorie neurali tra due condizioni sperimentali nello spazio delle 
% componenti principali (PC). L’obiettivo è valutare:
% - quanto le traiettorie dei due stati sperimentali differiscano in forma
%   per ciascun target comportamentale;
% - quanto la condizione 2 possa essere trasformata (rotazione, traslazione,
%   scala) per sovrapporsi alla condizione 1;
% - la similarità punto-per-punto delle traiettorie dopo l’allineamento.
%
% La procedura prevede:
% 1) normalizzazione Z-score delle due condizioni tramite media e deviazione 
%    standard globali (muZ, sigmaZ);
% 2) proiezione dei dati sulle componenti principali comuni (coeff);
% 3) selezione delle prime PC (pc_idx) per l’analisi delle traiettorie;
% 4) esecuzione dell’allineamento Procrustes target-per-target tra 
%    condizione 1 e condizione 2;
% 5) estrazione della distanza Procrustes d_t per ciascun target 
%    (d piccolo = forme simili);
% 6) conversione di d_t in uno score di similarità compreso tra 0 e 1;
% 7) visualizzazione:
%       - grafico a barre di d_t e dello score di similarità (1 - d_t),
%       - traiettorie 3D nelle prime PC prima/dopo l’allineamento,
%       - residuals temporali (distanza punto-per-punto) dopo Procrustes.
%
% Prima dell’esecuzione, assicurarsi di impostare correttamente:
%   condition_sep : cell array con le due condizioni da confrontare 
%                   (run codice pca_w_fixed.m)
%   muZ, sigmaZ   : parametri di normalizzazione globali 
%                   (run codice pca_w_fixed.m)
%   coeff         : componenti principali comuni ottenute da PCA condivisa 
%                   (run codice pca_w_fixed.m)
%   pc_idx        : indici delle PC da utilizzare per Procrustes
%   n_targets     : numero di target comportamentali
%   nbin          : numero di bin temporali per target
%   w             : finestra di smoothing per smoothdata
% =========================================================

%% Procrustes alignment
% Allineamento Procrustes tra le due condizioni.
% Per ogni target si trova la trasformazione (rotazione, traslazione e scala)
% che rende la traiettoria della condizione 2 il più simile possibile
% a quella della condizione 1. Un valore d piccolo indica forme simili.

Xz_cond1 = (condition_sep{1} - muZ) ./ sigmaZ;
Xz_cond2 = (condition_sep{2} - muZ) ./ sigmaZ;

scores1 = Xz_cond1 * coeff;   
scores2 = Xz_cond2 * coeff;

% Parametri per il confronto Procrustes
pc_idx = 1:3;                     % PC da usare 
T = nbin;   

% Seleziona solo le PC desiderate
scores1_pc = scores1(:, pc_idx);     
scores2_pc = scores2(:, pc_idx);     

% Procrustes per ogni target
d_targets = zeros(n_targets,1);
for tg = 1:n_targets
    % righe di X/Y corrispondenti al target tg
    idx = (tg-1)*T + (1:T);    

    Xt = scores1_pc(idx, :);     
    Yt = scores2_pc(idx, :);      

    % Procrustes per il singolo target
    [d_t, Yt_aligned, transform_t] = procrustes(Xt, Yt, ...
        'scaling', true, 'reflection', false);

    d_targets(tg) = d_t;
    traj_cond1{tg}  = Xt;
    traj_cond2A{tg} = Yt_aligned;
end
disp('Distanza di Procrustes per target:');
disp(d_targets);

%% Figure (1)
% Grafico a barre del valore d per ciascun target.
% Più è basso, più le due traiettorie hanno forma simile.

figure('Color', 'w')
bar(d_targets, 'FaceColor', [0.6 0.8 0.6], 'EdgeColor', 'none'); xlabel('Target');
ylim([0, max(1, max(d_targets))])
ylabel('Procrustes distance');
title('Similarity of trajectory shape between conditions');
box("off")


%% Figure (2)
% Conversione della distanza d in uno score di similarità compreso tra 0 e 1,
% utile per visualizzare quanto le traiettorie sono simili (1 = identiche).

similarity = 1 - d_targets;

figure('Color', 'w'); 
bar(similarity, 'FaceColor', [0.6 0.8 0.6], 'EdgeColor', 'none'); xlabel('Target');
ylim([0 1]);
ylabel('Similarity (1 - d)');
title('Similarity score between conditions');
box("off")


%% Figure (3) - Traiettorie dopo Procrustes
% Visualizzazione delle traiettorie nelle prime 3 PC.
% Per ogni target si mostra:
% - la traiettoria della condizione 1
% - la traiettoria della condizione 2 dopo Procrustes
% Le due curve dovrebbero sovrapporsi bene se l’allineamento è efficace.

figure('Color','w'); 
hold on;
axis equal; 
grid on;
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('Trajectories after Procrustes (Cond1 vs aligned Cond2)');

for tg = 1:n_targets
    Xt  = traj_cond1{tg};
    Ya  = traj_cond2A{tg};

    t1 = smoothdata(Xt(:,1), 'gaussian', w);
    t2 = smoothdata(Xt(:,2), 'gaussian', w);
    t3 = smoothdata(Xt(:,3), 'gaussian', w);

    a1 = smoothdata(Ya(:,1), 'gaussian', w);
    a2 = smoothdata(Ya(:,2), 'gaussian', w);
    a3 = smoothdata(Ya(:,3), 'gaussian', w);


    % traiettoria condizione 1 (linea continua)
    plot3(t1, t2, t3, '-', 'Color', colors(tg,:), 'LineWidth', 2, 'DisplayName', sprintf('Target %d - cond1', tg));

    % traiettoria condizione 2 allineata (linea tratteggiata)
    plot3(a1, a2, a3, '--', 'Color', colors(tg,:), 'LineWidth', 2, 'DisplayName', sprintf('Target %d - cond2 aligned', tg));

    % marker inizio/fine per condizione 1
    plot3(t1(1), t2(1), t3(1), 'o', ...
              'MarkerSize', 8, 'MarkerFaceColor', colors(tg,:), ...
              'MarkerEdgeColor', colors(tg,:), 'HandleVisibility','off');
    plot3(t1(end), t2(end), t3(end), '^', ...
              'MarkerSize', 8, 'MarkerFaceColor', colors(tg,:), ...
              'MarkerEdgeColor', colors(tg,:), 'HandleVisibility','off');

    % marker inizio/fine per condizione 2 allineata
    plot3(a1(1), a2(1), a3(1), 'o', ...
              'MarkerSize', 8, 'MarkerFaceColor', colors(tg,:), ...
              'MarkerEdgeColor', colors(tg,:), 'HandleVisibility','off');
    plot3(a1(end), a2(end), a3(end), '^', ...
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

legend([hStyles hStart hEnd], {condNames{:}, 'Start', 'End'}, ...
       'Location','northeastoutside');

%% Figure (4) - Residuals per ciascun target (quanto differiscono nel tempo)
% Per ogni target si calcola la distanza punto-per-punto tra le due
% traiettorie dopo Procrustes. Il plot sovrappone i residuals di tutti i target
% per mostrare in quali momenti le condizioni differiscono maggiormente.

figure('Color','w'); 
hold on; 
grid on;

for tg = 1:n_targets
    Xt = traj_cond1{tg};
    Ya = traj_cond2A{tg};
    residuals = sqrt(sum((Xt - Ya).^2, 2));

    plot(residuals, ...
         'LineWidth', 1.5, ...
         'Color', colors(tg,:), 'DisplayName', sprintf('%d', tg));
end

xlabel('Time bin');
ylabel('Residual distance');
title('Residuals after Procrustes (all targets)');
legend('Location','northeastoutside');