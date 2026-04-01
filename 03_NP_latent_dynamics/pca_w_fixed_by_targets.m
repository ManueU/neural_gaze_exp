% =========================================================
% DESCRIZIONE:
% Esegue un’analisi PCA su più condizioni sperimentali per valutare:
% - la separabilità dei target di movimento nello spazio delle componenti principali (PC);
% - la relazione tra le diverse condizioni nello stesso spazio PC.
%
% La finestra temporale analizzata è fissa rispetto all’onset del reach:
% w = [-100, +500] ms.
%
% Prima dell’esecuzione, modificare i parametri:
%   filename : nome del file .mat da caricare
%   PRE, POST: etichette delle condizioni di interesse
% =========================================================

% 1: medial arm 
% 2: lateral hand 
 
clearvars -except mean_baseline_common std_baseline_common 
close all
clc

sets = [1,2,3,4,5,6]; 
n_sets = numel(sets); 
n_arrays = 2;      
n_channels = 96;    
n_targets = 8; 
n_trials = 16; 
bin_size = 0.02; 
period_pre = 1.0; 
period_post = 0.5; 

filename = {'../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat',...
            '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat',...
            '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat',...
            '../00_Data_extraction/BCI02_Session_0924/gaze_BCI02_exclUpdated.mat'};


%% Costruzione matrice PCA 
rng(0)
condition = []; 
for d = 1:numel(filename) 
    fprintf('\nDataset: %s\n', filename{d}); 
    load(filename{d});

    [~, baseName, ext] = fileparts(filename{d});
    ds_name = [baseName ext];

    if strcmp(ds_name, 'controlled_BCI02_exclUpdated.mat')
        PRE = "Gaze";
        POST = "Reach";
    elseif strcmp(ds_name, 'gaze_BCI02_exclUpdated.mat')
        PRE = "Pres12";
        POST = "Gaze";
    else
        PRE = "Pres12";
        POST = "Reach";
    end
        
    idx_pres = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == PRE); 
    idx_reach = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == POST); 
    if isempty(idx_pres) || isempty(idx_reach)
        error('Stati PRE/POST non trovati: controlla PRE/POST');
    end

    pre_bins  = max(1, round(period_pre/bin_size));
    post_bins = max(1, round(period_post/bin_size));

    pca_matrix_array = []; 
    for array = 1:n_arrays
        pca_matrix = []; 
    
        for channel = 1:n_channels
            zscored_by_targets = []; 
            ch_global = (array-1)*n_channels + channel;

            for target = 1:n_targets
                firing_rate = []; 
                for set_ = 1:n_sets
                    set = sets(set_);

                    allT  = [data(set).Data(array).Interp.Target_ID];
                    allEx = [data(set).Data(array).Interp.Excluded];
                    
                    idx = find(allT == target & allEx == 0);

                    % nUse = min(15, numel(idx));  
                    % disp(nUse)
                    % idx  = idx(randperm(numel(idx), nUse));
                    for j = 1:length(idx) 
                        tmp_pre = data(set).Data(array).Interp(idx(j)).Task_states{idx_pres, 2}(end-pre_bins+1:end, channel); 
                        tmp_post = data(set).Data(array).Interp(idx(j)).Task_states{idx_reach,2}(1:post_bins, channel); 
                        vect = [tmp_pre; tmp_post];  
                        firing_rate = [firing_rate, vect./bin_size];  
                    end
                end
                zscored = (mean(firing_rate,2) - mean_baseline_common(ch_global))./std_baseline_common(ch_global);
                zscored(isnan(zscored) | isinf(zscored)) = 0;
                % zMax = 5;   % oppure 4 o 6
                % zscored(zscored >  zMax) =  zMax;
                % zscored(zscored < -zMax) = -zMax;
                zscored_by_targets = [zscored_by_targets; zscored]; 
                
            end 
            pca_matrix = [pca_matrix, zscored_by_targets]; 
        end
        pca_matrix_array = [pca_matrix_array, pca_matrix]; 
    end 
    condition = [condition; pca_matrix_array];
    condition_sep(d) = {pca_matrix_array};
end 


%% PCA
[coeff, score, latent, tsquared, explained] = pca(condition, 'Algorithm','svd');

%% Varianza spiegata
fprintf('\nVarianza spiegata:\n');
for k = [1 3 5 10 20 50]
    if k <= numel(explained)
        fprintf('  Prime %2d PC: %5.2f %%\n', k, sum(explained(1:k)));
    end
end
fprintf('  Totale     : %5.2f %%\n\n', sum(explained));

figure('Color','w');
subplot(2,1,1)
bar(explained)
xlim([1 size(score,2)])
xlabel('Principal Component')
ylabel('Explained Variance (%)')
title('Explained Variance per Principal Component')
grid on

% Varianza spiegata cumulativa
subplot(2,1,2)
plot(cumsum(explained), '-o','LineWidth',1.5)
xlabel('Number of Components')
ylabel('Cumulative Explained Variance (%)')
title('Cumulative Explained Variance')
xlim([1 size(score,2)])
ylim([0 100])
grid on


%% Parameters for figures
styles = {'-','--',':','-.'};
colors = [
    0.839, 0.153, 0.157;  % rosso
    0.122, 0.467, 0.706;  % blu
    0.172, 0.627, 0.172;  % verde
    0.580, 0.404, 0.741;  % viola
    1.000, 0.498, 0.055;  % arancione
    0.737, 0.741, 0.133;  % giallo oliva
    0.549, 0.337, 0.294;  % marrone
    0.890, 0.466, 0.760;  % rosa
];

nCond = numel(filename);
T = zeros(nCond,1);
for c=1:nCond
    T(c) = size(condition_sep{c},1);
end
offset = [0; cumsum(T(1:end-1))];
labels = {'Free-gaze', 'Gaze-on-center', 'Gaze-on-target', 'Gaze-only'};
w = 11; 

%% Figure (1) - plot dell'intero intervallo
figure('Color','White'); 
hold on

for c = 1:nCond
    nbin = size(condition_sep{c},1) / n_targets; 
    for t = 1:n_targets
        rows = offset(c) + (t-1)*nbin + (1:nbin);
        traj = score(rows,1:3);
        t1 = smoothdata(traj(:,1), 'gaussian', w);
        t2 = smoothdata(traj(:,2), 'gaussian', w);
        t3 = smoothdata(traj(:,3), 'gaussian', w);

        plot3(t1, t2, t3, 'Color', colors(t,:), 'LineWidth', 2, 'LineStyle',styles{c});
        hold on
    
        % marker inizio (cerchio)
        plot3(t1(1), t2(1), t3(1), 'o', ...
              'MarkerSize', 8, 'MarkerFaceColor', colors(t,:), ...
              'MarkerEdgeColor', colors(t,:), 'HandleVisibility','off');
    
        % marker fine (triangolo)
        plot3(t1(end), t2(end), t3(end), '^', ...
              'MarkerSize', 8, 'MarkerFaceColor', colors(t,:), ...
              'MarkerEdgeColor', colors(t,:), 'HandleVisibility','off');
    end
end 
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
grid on; axis equal;

hStart = plot3(nan,nan,nan, 'o', 'Color','k', 'MarkerFaceColor','k');
hEnd   = plot3(nan,nan,nan, '^', 'Color','k', 'MarkerFaceColor','k');
if nCond == 1
    legend([hStart hEnd], {'Start','End'}, 'Location','northeastoutside');
else
    styles = {'-','--',':','-.'};
    hStyles = gobjects(1, nCond);
    for c = 1:nCond
        % Dummy per linea condizione
        ls = styles{mod(c-1, numel(styles))+1};
        hStyles(c) = plot3(nan,nan,nan, ls, 'Color','k', 'LineWidth', 2);
    end
    legEntries{nCond+1} = 'Start';
    legEntries{nCond+2} = 'End';
    legend([hStyles hStart hEnd], labels, 'Location','northeastoutside');
end


%% Figure (2) - plot dell'intervallo di interesse
bin_pre_keep  = round(0.1 / bin_size);  % 100 ms
bin_post_keep = post_bins;              % 500 ms
idx_pre  = (pre_bins - bin_pre_keep + 1) : pre_bins;
idx_post = pre_bins + (1:bin_post_keep);
idx_keep = [idx_pre idx_post];

figure('Color','White'); 
hold on

for c = 1:nCond
    nbin = size(condition_sep{c},1) / n_targets; 
    for t = 1:n_targets
        rows = offset(c) + (t-1)*nbin + (1:nbin);
        traj = score(rows,1:3);
        t1 = smoothdata(traj(:,1), 'gaussian', w);
        t2 = smoothdata(traj(:,2), 'gaussian', w);
        t3 = smoothdata(traj(:,3), 'gaussian', w);

        plot3(t1(idx_keep), t2(idx_keep), t3(idx_keep), 'Color', colors(t,:), 'LineWidth', 2, 'LineStyle', styles{c});
        hold on
    
        % start (-100 ms)
        plot3(t1(idx_keep(1)), t2(idx_keep(1)), t3(idx_keep(1)), 'o', ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', colors(t,:), ...
            'MarkerEdgeColor', colors(t,:), ...
            'HandleVisibility','off');
        
        % end (+500 ms)
        plot3(t1(idx_keep(end)), t2(idx_keep(end)), t3(idx_keep(end)), '^', ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', colors(t,:), ...
            'MarkerEdgeColor', colors(t,:), ...
            'HandleVisibility','off');
    end
end 
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
grid on; axis equal;

hStart = plot3(nan,nan,nan, 'o', 'Color','k', 'MarkerFaceColor','k');
hEnd   = plot3(nan,nan,nan, '^', 'Color','k', 'MarkerFaceColor','k');
if nCond == 1
    legend([hStart hEnd], {'Start','End'}, 'Location','northeastoutside');
else
    styles = {'-','--',':','-.'};
    hStyles = gobjects(1, nCond);
    for c = 1:nCond
        % Dummy per linea condizione
        ls = styles{mod(c-1, numel(styles))+1};
        hStyles(c) = plot3(nan,nan,nan, ls, 'Color','k', 'LineWidth', 2);
    end
    legEntries{nCond+1} = 'Start';
    legEntries{nCond+2} = 'End';
    legend([hStyles hStart hEnd], labels, 'Location','northeastoutside');
end

%% Figure (3) - PCA centroids (1 dot per target x condition)
bin_pre_keep  = round(0.1 / bin_size);  % 100 ms
bin_post_keep = post_bins;              % 500 ms
idx_pre  = (pre_bins - bin_pre_keep + 1) : pre_bins;
idx_post = pre_bins + (1:bin_post_keep);
idx_keep = [idx_pre idx_post];

figure('Color','white'); 
hold on
styles = {'o','s','^','d','v','>','<','p'};  
centroids = nan(nCond, n_targets, 3);

hCond = gobjects(1, nCond);
for c = 1:nCond
    nbin = size(condition_sep{c},1) / n_targets;
    mk = styles{mod(c-1, numel(styles)) + 1};

    for t = 1:n_targets
        rows = offset(c) + (t-1)*nbin + (1:nbin);
        traj = score(rows, 1:3);
        seg  = traj(idx_keep, :);  % solo -100..+500 ms
        
        mu = mean(seg, 1, 'omitnan');
        centroids(c,t,:) = mu;

        h = plot3(mu(1), mu(2), mu(3), mk, ...
            'MarkerSize', 10, ...
            'MarkerFaceColor', colors(t,:), ...
            'MarkerEdgeColor', 'none', ...
            'LineWidth', 0.8);

        % dummy handle (nero)
        hCond(c) = plot3(nan, nan, nan, mk, ...
            'MarkerSize', 10, ...
            'MarkerFaceColor', 'k', ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.8);
    end
end

legend(hCond, labels, 'Location','northeastoutside');

xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
grid on; axis equal; view(3);

%% Figure (4) === PCA centroids: emphasize condition separation on PC1-PC2 ===

figure('Color','white'); 
hold on; grid on; box on

condColor = [
    0.6  0.6  0.6;   % free-gaze   (gray)
    0.05 0.10 0.35;  % motor       (navy blue)
    0.55 0.10 0.20;  % controlled  (bordeaux)
    0 0 0;           % gaze-only  (black)
];

% Marker per TARGET (uno per target) - se vuoi puoi anche usare text label
targetMarkers = {'o','s','d','^','v','>','<','p'};  % per 8 target
ms = 10;

% raccogli punti PC1-PC2
X = nan(nCond*n_targets,2);
G = strings(nCond*n_targets,1); % label condizione
idx = 0;

condNames = strings(nCond,1);
for c = 1:nCond
    [~, baseName, ~] = fileparts(filename{c});
    [namePart, ~] = strtok(baseName, '_');
    condNames(c) = string(namePart);
end

for c = 1:nCond
    for t = 1:n_targets
        idx = idx + 1;
        mu = squeeze(centroids(c,t,1:2))';  % PC1-PC2
        X(idx,:) = mu;
        G(idx) = condNames(c);

        mk = targetMarkers{t};
        plot(mu(1), mu(2), mk, ...
            'MarkerSize', ms, ...
            'MarkerFaceColor', condColor(c,:), ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.8, ...
            'HandleVisibility','off');
    end

    % === centroide condizione (media dei target) ===
    muC = mean(squeeze(centroids(c,:,1:2)), 1, 'omitnan');

    plot(muC(1), muC(2), 'x', ...
        'Color', condColor(c,:), ...
        'LineWidth', 2.5, ...
        'MarkerSize', 12);

    % === confidence ellipse (1 SD o 95%) ===
    pts = squeeze(centroids(c,:,1:2));          % [targets x 2]
    pts = pts(all(isfinite(pts),2),:);
    if size(pts,1) >= 3
        C = cov(pts);
        [V,D] = eig(C);
        % scala: 95% approx chi2inv(0.95,2)=5.991 (senza Stats TB: valore hardcoded)
        k = sqrt(5.991);  % 95%
        th = linspace(0,2*pi,200);
        ell = (V*sqrt(D))*[cos(th); sin(th)] * k;
        ell = ell' + muC;

        plot(ell(:,1), ell(:,2), '-', 'Color', condColor(c,:), 'LineWidth', 2);
    end
end

% legenda: una entry per condizione (dummy)
hLeg = gobjects(1,nCond);
for c = 1:nCond
    hLeg(c) = plot(nan,nan,'o', ...
        'MarkerSize', ms, ...
        'MarkerFaceColor', condColor(c,:), ...
        'MarkerEdgeColor', 'k');
end
legend(hLeg, labels, 'Location','northeastoutside');

xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
axis equal


%% Video
% outputVideo = VideoWriter('pca_gaze_BCI02.mp4', 'Uncompressed AVI'); 
% outputVideo.FrameRate = 15;  
% open(outputVideo);
% 
% for angle = 0:2:360
%     view(angle, 30);   
%     frame = getframe(gcf);  
%     writeVideo(outputVideo, frame);  
% end
% close(outputVideo);






