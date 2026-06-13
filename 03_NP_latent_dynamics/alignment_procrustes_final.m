%% =========================================================
% PROCRUSTES ALIGNMENT - ALL PAIRWISE COMPARISONS
% =========================================================

pc_idx = 1:3;
T = size(condition_sep{1},1) / n_targets;

nCond = numel(condition_sep);
cond_pairs = nchoosek(1:nCond, 2);
nPairs = size(cond_pairs,1);

condNames = { ...
    'Free-gaze', ...
    'Gaze-on-center', ...
    'Gaze-on-target' ...
    };

pair_names = cell(nPairs,1);
for p = 1:nPairs
    pair_names{p} = sprintf('%s vs %s', ...
        condNames{cond_pairs(p,1)}, ...
        condNames{cond_pairs(p,2)});
end

% Proiezione PCA di tutte le condizioni
scores_pc = cell(1,nCond);
for c = 1:nCond
    scores = (condition_sep{c} - mu_pca) * coeff;
    scores_pc{c} = scores(:, pc_idx);
end

% Output
D_targets = nan(nPairs, n_targets);
D_global  = nan(nPairs,1);
Residuals = cell(nPairs, n_targets);
TrajRef   = cell(nPairs, n_targets);
TrajAlign = cell(nPairs, n_targets);

%% =========================================================
% LOOP AUTOMATICO SU TUTTE LE COPPIE
% =========================================================

for p = 1:nPairs

    c1 = cond_pairs(p,1);
    c2 = cond_pairs(p,2);

    X_all = scores_pc{c1};
    Y_all = scores_pc{c2};

    % -------- Procrustes globale --------
    [D_global(p), Y_all_aligned] = procrustes(X_all, Y_all, ...
        'scaling', false, ...
        'reflection', false);

    % -------- Procrustes target-by-target --------
    for tg = 1:n_targets

        idx = (tg-1)*T + (1:T);

        Xt = X_all(idx,:);
        Yt = Y_all(idx,:);

        [d_t, Yt_aligned] = procrustes(Xt, Yt, ...
            'scaling', false, ...
            'reflection', false);

        D_targets(p,tg) = d_t;
        TrajRef{p,tg}   = Xt;
        TrajAlign{p,tg} = Yt_aligned;

        Residuals{p,tg} = sqrt(sum((Xt - Yt_aligned).^2, 2));
    end
end

%% =========================================================
% FIGURE 1 - DISTANZE TARGET-BY-TARGET PER OGNI COPPIA
% =========================================================

figure('Color','w');

for p = 1:nPairs
    subplot(1,nPairs,p)

    bar(D_targets(p,:), ...
        'FaceColor', [0.6 0.8 0.6], ...
        'EdgeColor', 'none');

    xlabel('Target');
    ylabel('Procrustes distance');
    title(pair_names{p}, 'Interpreter','none');

    ylim([0 max(1, max(D_targets(:)))]);
    box off
    grid on
end

sgtitle('Target-by-target Procrustes distance');

%% =========================================================
% FIGURE 2 - INDICE GLOBALE PAIRWISE
% =========================================================
clear set
figure('Color','w');

bar(D_global, ...
    'FaceColor', [0.6 0.8 0.6], ...
    'EdgeColor', 'none');

set(gca, ...
    'XTick', 1:nPairs, ...
    'XTickLabel', pair_names, ...
    'XTickLabelRotation', 30);

ylabel('Global Procrustes distance');
title('Global trajectory similarity between conditions');
box off
grid on

%% =========================================================
% FIGURE 3 - RESIDUALS TEMPORALI PER OGNI COPPIA
% =========================================================

figure('Color','w');

for p = 1:nPairs
    subplot(1,nPairs,p)
    hold on
    grid on

    for tg = 1:n_targets
        plot(Residuals{p,tg}, ...
            'LineWidth', 1.5, ...
            'Color', colors(tg,:), ...
            'DisplayName', sprintf('T%d', tg));
    end

    xlabel('Time bin');
    ylabel('Residual distance');
    title(pair_names{p}, 'Interpreter','none');
    box off
end

sgtitle('Residuals after Procrustes');

%% =========================================================
% FIGURE 4 - TRAIETTORIE ALLINEATE PER UNA COPPIA SELEZIONATA
% =========================================================

pair_to_plot = 1;   % cambia 1, 2, 3 per scegliere la coppia

figure('Color','w');
hold on
axis equal
grid on

xlabel('PC1');
ylabel('PC2');
zlabel('PC3');

title(sprintf('Trajectories after Procrustes: %s', ...
    pair_names{pair_to_plot}), ...
    'Interpreter','none');

for tg = 1:n_targets

    Xt = TrajRef{pair_to_plot,tg};
    Ya = TrajAlign{pair_to_plot,tg};

    t1 = smoothdata(Xt(:,1), 'gaussian', w);
    t2 = smoothdata(Xt(:,2), 'gaussian', w);
    t3 = smoothdata(Xt(:,3), 'gaussian', w);

    a1 = smoothdata(Ya(:,1), 'gaussian', w);
    a2 = smoothdata(Ya(:,2), 'gaussian', w);
    a3 = smoothdata(Ya(:,3), 'gaussian', w);

    plot3(t1, t2, t3, '-', ...
        'Color', colors(tg,:), ...
        'LineWidth', 2);

    plot3(a1, a2, a3, '--', ...
        'Color', colors(tg,:), ...
        'LineWidth', 2);

    plot3(t1(1), t2(1), t3(1), 'o', ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', colors(tg,:), ...
        'MarkerEdgeColor', colors(tg,:), ...
        'HandleVisibility','off');

    plot3(t1(end), t2(end), t3(end), '^', ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', colors(tg,:), ...
        'MarkerEdgeColor', colors(tg,:), ...
        'HandleVisibility','off');
end

h1 = plot3(nan,nan,nan,'-','Color','k','LineWidth',2);
h2 = plot3(nan,nan,nan,'--','Color','k','LineWidth',2);

legend([h1 h2], ...
    {condNames{cond_pairs(pair_to_plot,1)}, ...
     [condNames{cond_pairs(pair_to_plot,2)} ' aligned']}, ...
    'Location','northeastoutside');

%% =========================================================
% FIGURE - RESIDUALS MEDIATI SUI TARGET ± SEM
% =========================================================

figure('Color','w');
hold on
grid on
box off

time_ms = linspace(-100, 500, T);
x = time_ms;

pair_colors = [
    0.2 0.2 0.2;
    0.4 0.4 0.8;
    0.8 0.4 0.4
];


for p = 1:nPairs

    R = nan(n_targets, T);

    for tg = 1:n_targets
        R(tg,:) = Residuals{p,tg};
    end

    muR  = mean(R, 1, 'omitnan');
    semR = std(R, [], 1, 'omitnan') ./ sqrt(sum(~isnan(R),1));

    % Area SEM
    fill([x fliplr(x)], ...
         [muR + semR, fliplr(muR - semR)], ...
         pair_colors(p,:), ...
         'FaceAlpha', 0.18, ...
         'EdgeColor', 'none', ...
         'HandleVisibility','off');
    
    
    % Media
    plot(x, muR, ...
        'Color', pair_colors(p,:), ...
        'LineWidth', 2.5, ...
        'DisplayName', pair_names{p});
end

xlabel('Time from reach onset (ms)');
ylabel('Residual distance');
title('Mean residuals after Procrustes across targets');
legend('Location','northeastoutside');