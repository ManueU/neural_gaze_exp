%% =========================================================
% CCA AUTOMATICA SU TUTTE LE COPPIE DI CONDIZIONI
% Richiede già disponibili:
% condition_sep, coeff, mu_pca, filename, n_targets, colors, w
%% =========================================================

nCond = numel(condition_sep);
condPairs = nchoosek(1:nCond, 2);

m_values = [4 6 8 10 12 15 20];
m_final  = 10;
nCC = 4;
nBoot = 20;
r_thr = 0.5;

styles = {'-','--'};
condNames = string(labels(1:nCond));

CCA_results = struct();

for ipair = 1:size(condPairs,1)

    c1 = condPairs(ipair,1);
    c2 = condPairs(ipair,2);

    fprintf('\n============================================\n');
    fprintf('CCA: %s vs %s\n', condNames(c1), condNames(c2));
    fprintf('============================================\n');

    %% Normalizzazione rispetto alla PCA comune
    Xz_cond1 = condition_sep{c1} - mu_pca;
    Xz_cond2 = condition_sep{c2} - mu_pca;

    %% Proiezione sulle PC comuni
    scores1 = Xz_cond1 * coeff;
    scores2 = Xz_cond2 * coeff;

    %% Controllo dimensioni
    if size(scores1,1) ~= size(scores2,1)
        error('Le condizioni %d e %d hanno numero di righe diverso.', c1, c2);
    end

    nbin = size(scores1,1) / n_targets;

    if mod(size(scores1,1), n_targets) ~= 0
        error('Il numero di righe non è divisibile per n_targets.');
    end

    %% =====================================================
    % CCA stability
    %% =====================================================

    r_all = nan(numel(m_values), nCC, nBoot);

    for im = 1:numel(m_values)

        m = m_values(im);
        pc_idx_cca = 1:m;

        Xcca_full = scores1(:, pc_idx_cca);
        Ycca_full = scores2(:, pc_idx_cca);

        nSamples = size(Xcca_full, 1);

        for b = 1:nBoot

            idx = randsample(nSamples, nSamples, true);

            Xcca = zscore(Xcca_full(idx, :));
            Ycca = zscore(Ycca_full(idx, :));

            [~, ~, r, ~, ~, ~] = canoncorr(Xcca, Ycca);

            nSave = min(nCC, numel(r));
            r_all(im, 1:nSave, b) = r(1:nSave);

        end
    end

    r_mean = mean(r_all, 3, 'omitnan');
    r_std  = std(r_all, 0, 3, 'omitnan');

    fprintf('Canonical correlations, prime %d CC:\n', nCC);

    for im = 1:numel(m_values)
        fprintf('m = %2d:  ', m_values(im));
        for k = 1:nCC
            fprintf('CC%d = %.3f ± %.3f   ', k, r_mean(im,k), r_std(im,k));
        end
        fprintf('\n');
    end

    %% Figure 1 - stability
    figure('Color','w');
    hold on;

    for k = 1:nCC
        errorbar(m_values, r_mean(:,k), r_std(:,k), '-o', ...
            'DisplayName', sprintf('CC%d', k));
    end

    ylim([0 1]);
    xlabel('PCs');
    ylabel('Canonical correlation');
    title(sprintf('CCA stability: %s vs %s', condNames(c1), condNames(c2)));
    legend('Location','best');
    grid on;

    %% =====================================================
    % CCA finale con m_final
    %% =====================================================

    pc_idx_cca = 1:m_final;

    Xcca = zscore(scores1(:, pc_idx_cca));
    Ycca = zscore(scores2(:, pc_idx_cca));

    [A, B, r, U, V, stats] = canoncorr(Xcca, Ycca);

    disp('Canonical correlations finali:');
    disp(r);

    k_F = min(nCC, numel(r));
    d_CCA_final = 1 - mean(r(1:k_F));

    fprintf('CCA-distance, prime %d CC, m = %d: %.4f\n', ...
        k_F, m_final, d_CCA_final);

    %% Salvataggio risultati
    CCA_results(ipair).pair = [c1 c2];
    CCA_results(ipair).condNames = [condNames(c1), condNames(c2)];
    CCA_results(ipair).m_values = m_values;
    CCA_results(ipair).r_all = r_all;
    CCA_results(ipair).r_mean = r_mean;
    CCA_results(ipair).r_std = r_std;
    CCA_results(ipair).A = A;
    CCA_results(ipair).B = B;
    CCA_results(ipair).r = r;
    CCA_results(ipair).U = U;
    CCA_results(ipair).V = V;
    CCA_results(ipair).stats = stats;
    CCA_results(ipair).d_CCA_final = d_CCA_final;

    %% =====================================================
    % Figure 2 - traiettorie nello spazio canonico
    %% =====================================================

    figure('Color','w');
    hold on;
    axis equal;
    grid on;
    view(3);

    xlabel('Can1');
    ylabel('Can2');
    zlabel('Can3');

    title(sprintf('CCA trajectories: %s vs %s', ...
        condNames(c1), condNames(c2)));

    for tg = 1:n_targets

        idx = (tg-1)*nbin + (1:nbin);

        Ut = U(idx, 1:3);
        Vt = V(idx, 1:3);

        u1 = smoothdata(Ut(:,1), 'gaussian', w);
        u2 = smoothdata(Ut(:,2), 'gaussian', w);
        u3 = smoothdata(Ut(:,3), 'gaussian', w);

        v1 = smoothdata(Vt(:,1), 'gaussian', w);
        v2 = smoothdata(Vt(:,2), 'gaussian', w);
        v3 = smoothdata(Vt(:,3), 'gaussian', w);

        plot3(u1, u2, u3, '-', ...
            'Color', colors(tg,:), ...
            'LineWidth', 2, ...
            'HandleVisibility','off');

        plot3(v1, v2, v3, '--', ...
            'Color', colors(tg,:), ...
            'LineWidth', 2, ...
            'HandleVisibility','off');

        plot3(u1(1), u2(1), u3(1), 'o', ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', colors(tg,:), ...
            'MarkerEdgeColor', colors(tg,:), ...
            'HandleVisibility','off');

        plot3(u1(end), u2(end), u3(end), '^', ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', colors(tg,:), ...
            'MarkerEdgeColor', colors(tg,:), ...
            'HandleVisibility','off');

        plot3(v1(1), v2(1), v3(1), 'o', ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', colors(tg,:), ...
            'MarkerEdgeColor', colors(tg,:), ...
            'HandleVisibility','off');

        plot3(v1(end), v2(end), v3(end), '^', ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', colors(tg,:), ...
            'MarkerEdgeColor', colors(tg,:), ...
            'HandleVisibility','off');
    end

    h1 = plot3(nan,nan,nan, '-',  'Color','k', 'LineWidth', 2);
    h2 = plot3(nan,nan,nan, '--', 'Color','k', 'LineWidth', 2);
    hStart = plot3(nan,nan,nan, 'o', 'Color','k', 'MarkerFaceColor','k');
    hEnd   = plot3(nan,nan,nan, '^', 'Color','k', 'MarkerFaceColor','k');

    legend([h1 h2 hStart hEnd], ...
        {char(condNames(c1)), char(condNames(c2)), 'Start', 'End'}, ...
        'Location','northeastoutside');

    %% =====================================================
    % Figure 3 - canonical correlations
    %% =====================================================

    figure('Color','w');
    bar(r, 'FaceColor', [0.6 0.8 0.6], 'EdgeColor','none');
    ylim([0 1]);
    xlabel('Canonical dimension');
    ylabel('Correlation');
    title(sprintf('Canonical correlations: %s vs %s', ...
        condNames(c1), condNames(c2)));
    box off;

    %% =====================================================
    % Figure 4 - canonical correlations con soglia
    %% =====================================================

    figure('Color','w');
    bar(r, 'FaceColor', [0.6 0.8 0.6], 'EdgeColor','none');
    yline(r_thr, '--', 'Shared threshold');
    ylim([0 1]);
    xlabel('Canonical dimension');
    ylabel('Correlation');
    title(sprintf('Shared dimensions: %s vs %s', ...
        condNames(c1), condNames(c2)));
    box off;

end