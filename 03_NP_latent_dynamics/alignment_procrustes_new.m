%% Procrustes alignment - trasformazione GLOBALE


Xz_cond1 = (condition_sep{1} - muZ);
% Xz_cond1 = (condition_sep{1} - muZ) ./ sigmaZ;

Xz_cond2 = (condition_sep{2} - muZ);
% Xz_cond2 = (condition_sep{2} - muZ) ./ sigmaZ;

% Proiezione sulle PC comuni
scores1 = Xz_cond1 * coeff;   
scores2 = Xz_cond2 * coeff;

% Parametri per il confronto Procrustes
pc_idx = 1:3;      % PC da usare 
T      = nbin;     % numero di bin per target

% Seleziona solo le PC desiderate
scores1_pc = scores1(:, pc_idx);     
scores2_pc = scores2(:, pc_idx);     

%% Procrustes complessivo (unico per tutte le traiettorie)

[d_global, Y_all_aligned, transform_all] = procrustes(scores1_pc, scores2_pc, ...
    'scaling', true, 'reflection', false);

disp('Distanza di Procrustes globale (tutti i target):');
disp(d_global);

similarity_global = 1 - d_global;
disp('Similarità globale (1 - d):');
disp(similarity_global);

%% Ricostruzione traiettorie per target usando la trasformazione GLOBALE

traj_cond1  = cell(n_targets,1);
traj_cond2A = cell(n_targets,1);

for tg = 1:n_targets
    % righe corrispondenti al target tg
    idx = (tg-1)*T + (1:T);
    
    % traiettorie nelle prime PC
    traj_cond1{tg}  = scores1_pc(idx, :);      % condizione 1
    traj_cond2A{tg} = Y_all_aligned(idx, :);   % condizione 2 ALLINEATA GLOBALMENTE
end

%% (Opzionale) misura per-target derivata dall'allineamento globale
% Qui NON si rifà Procrustes: si calcola una distanza media residua per target
d_targets_global = zeros(n_targets,1);
for tg = 1:n_targets
    Xt = traj_cond1{tg};
    Ya = traj_cond2A{tg};
    % distanza media normalizzata per target (puoi cambiare metrica se preferisci)
    d_targets_global(tg) = mean(sqrt(sum((Xt - Ya).^2, 2)));
end

disp('Distanza residua media per target (dopo Procrustes globale):');
disp(d_targets_global);

%% Figure (1) - esempio: bar plot delle distanze per target (globali)

figure('Color', 'w')
bar(d_targets_global, 'FaceColor', [0.6 0.8 0.6], 'EdgeColor', 'none'); 
xlabel('Target');
ylabel('Mean residual distance (global Procrustes)');
title('Residual per target dopo Procrustes GLOBALE');
box off

%% Figure (3) - Traiettorie dopo Procrustes GLOBALE

figure('Color','w'); 
hold on;
axis equal; 
grid on;
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title('Trajectories after GLOBAL Procrustes (Cond1 vs aligned Cond2)');

for tg = 1:n_targets
    Xt = traj_cond1{tg};
    Ya = traj_cond2A{tg};

    % smoothing
    t1 = smoothdata(Xt(:,1), 'gaussian', w);
    t2 = smoothdata(Xt(:,2), 'gaussian', w);
    t3 = smoothdata(Xt(:,3), 'gaussian', w);

    a1 = smoothdata(Ya(:,1), 'gaussian', w);
    a2 = smoothdata(Ya(:,2), 'gaussian', w);
    a3 = smoothdata(Ya(:,3), 'gaussian', w);

    % traiettoria condizione 1 (linea continua)
    plot3(t1, t2, t3, '-', 'Color', colors(tg,:), 'LineWidth', 2, ...
          'DisplayName', sprintf('Target %d - cond1', tg));

    % traiettoria condizione 2 allineata GLOBALMENTE (linea tratteggiata)
    plot3(a1, a2, a3, '--', 'Color', colors(tg,:), 'LineWidth', 2, ...
          'DisplayName', sprintf('Target %d - cond2 aligned (global)', tg));

    % marker inizio/fine per condizione 1
    plot3(t1(1), t2(1), t3(1), 'o', 'MarkerSize', 8, ...
          'MarkerFaceColor', colors(tg,:), 'MarkerEdgeColor', colors(tg,:), ...
          'HandleVisibility','off');
    plot3(t1(end), t2(end), t3(end), '^', 'MarkerSize', 8, ...
          'MarkerFaceColor', colors(tg,:), 'MarkerEdgeColor', colors(tg,:), ...
          'HandleVisibility','off');

    % marker inizio/fine per condizione 2 allineata
    plot3(a1(1), a2(1), a3(1), 'o', 'MarkerSize', 8, ...
          'MarkerFaceColor', colors(tg,:), 'MarkerEdgeColor', colors(tg,:), ...
          'HandleVisibility','off');
    plot3(a1(end), a2(end), a3(end), '^', 'MarkerSize', 8, ...
          'MarkerFaceColor', colors(tg,:), 'MarkerEdgeColor', colors(tg,:), ...
          'HandleVisibility','off');
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
% legend([hStyles hStart hEnd], {condNames{:}, 'Start', 'End'}, ...
%        'Location','northeastoutside');

%% Figure (4) - Residuals dopo Procrustes GLOBALE

figure('Color','w'); 
hold on; 
grid on;

for tg = 1:n_targets
    Xt = traj_cond1{tg};
    Ya = traj_cond2A{tg};
    residuals = sqrt(sum((Xt - Ya).^2, 2));

    plot(residuals, 'LineWidth', 1.5, ...
         'Color', colors(tg,:), 'DisplayName', sprintf('%d', tg));
end

xlabel('Time bin');
ylabel('Residual distance');
title('Residuals after GLOBAL Procrustes (all targets)');
legend('Location','northeastoutside');
