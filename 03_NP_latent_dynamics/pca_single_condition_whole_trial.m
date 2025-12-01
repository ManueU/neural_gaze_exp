clearvars 
close all
clc

sets_pca = [1,2,4,5,6];
n_sets = 5; 
n_arrays = 2;      
n_channels = 96;    
n_targets = 8; 
bin_size = 0.02; 


%% Load 
filename = '../00_Data_extraction/controlled_BCI02.mat';
load(filename)

%% Costruzione matrice PCA 
pca_matrix_array = []; 
for array = 1:n_arrays
    pca_matrix = []; 

    for channel = 1:n_channels
        firing_rate = []; 
        
        for target = 1:n_targets
            M_spikes = [];

            for set_all = 1:n_sets
                set = sets_pca(set_all);
                idx = find([data(set).Data(array).Interp.Target_ID] == target); 

                for j = 1:length(idx)
                    matrix = data(set).Data(array).Interp(idx(j)).Trial(1:460,channel); 
                    
                    M_spikes = [M_spikes, matrix];
                end

            end 
            M_spikes_mean = mean(M_spikes, 2);     
            firing_rate = [firing_rate; M_spikes_mean ./ bin_size];

        end 
        pca_matrix = [pca_matrix, firing_rate]; 

    end
    pca_matrix_array = [pca_matrix_array, pca_matrix]; 

end 

%% Filtraggio neuroni poco informativi
% varianza e firing rate medio per neurone
meanFR = mean(pca_matrix_array, 1);          % Hz
varFR  = var(pca_matrix_array, 0, 1);

minMeanFR = 1;       % Hz: neuroni con FR medio < 1 Hz vengono scartati
minVarFR  = 0.5;     % varianza minima (su firing rate in Hz)

valid_cols = (meanFR > minMeanFR) & (varFR > minVarFR);
fprintf('\nNeuroni totali: %d, neuroni tenuti dopo filtro: %d\n', numel(meanFR), sum(valid_cols));

pca_matrix_array = pca_matrix_array(:, valid_cols);

%% PCA
[Xz, muZ, sigmaZ] = zscore(pca_matrix_array, 0, 1); 
[coeff, score, latent, tsquared, explained] = pca(Xz, 'Algorithm','svd');

%% Explained variance
fprintf('\nVarianza spiegata:\n');
for k = [1 3 5 10 20 50]
    if k <= numel(explained)
        fprintf('  Prime %2d PC: %5.2f %%\n', k, sum(explained(1:k)));
    end
end
fprintf('  Totale     : %5.2f %%\n\n', sum(explained));

maxPCs = min(50, size(coeff,2));
cumExplained = cumsum(explained);
cumExplained = cumExplained(1:maxPCs);


%% Figure (1)
% figure('Color','white'); hold on
% 
% plot(1:maxPCs, cumExplained, 'b', 'LineWidth', 1);
% xlabel('PCs');
% ylabel('Cum. neural var. expl. (%)');
% xlim([1 maxPCs]);
% ylim([0 100]);
% grid on;
% box off;
% title('Cumulative neural variance explained');

%% Figure (2) - PCA
nbin = size(data(1).Data(1).Interp(1).Trial,1);   
nbin = 460; 
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

figure('Color','White')
w = 15; 
for t = 1:n_targets
    subplot(2,4,t)
    idx = (t-1)*nbin + (1:nbin);   
    traj = score(idx,1:3);        
    t1 = smoothdata(traj(:,1), 'gaussian', w);
    t2 = smoothdata(traj(:,2), 'gaussian', w);
    t3 = smoothdata(traj(:,3), 'gaussian', w);

    plot3(t1, t2, t3, ...
          'Color', colors(t,:), 'LineWidth', 2, 'MarkerFaceColor', colors(t,:));
    hold on

    % marker inizio (cerchio)
    plot3(t1(1), t2(1), t3(1), 'o', ...
          'MarkerSize', 8, 'MarkerFaceColor', colors(t,:), ...
          'MarkerEdgeColor', colors(t,:), 'HandleVisibility','off');
    hold on

    % marker fine (triangolo)
    plot3(t1(end), t2(end), t3(end), '^', ...
          'MarkerSize', 8, 'MarkerFaceColor', colors(t,:), ...
          'MarkerEdgeColor', colors(t,:), 'HandleVisibility','off');

    xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
    ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
    zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
    grid on; axis equal;

end

h1 = plot3(nan, nan, nan, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
h2 = plot3(nan, nan, nan, '^', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% legend([h1 h2], {'Start', 'End'}, 'Location', 'northeastoutside');


%% Figure (3) - PCA
nbin = size(data(1).Data(1).Interp(1).Trial,1);     
nbin = 460; 

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

figure('Color','White')
w = 40; 
for t = 1:n_targets
    idx = (t-1)*nbin + (1:nbin);   
    traj = score(idx,1:3);        
    t1 = smoothdata(traj(:,1), 'gaussian', w);
    t2 = smoothdata(traj(:,2), 'gaussian', w);
    t3 = smoothdata(traj(:,3), 'gaussian', w);

    plot3(t1, t2, t3, ...
          'Color', colors(t,:), 'LineWidth', 2, 'MarkerFaceColor', colors(t,:));
    hold on

    % marker inizio (cerchio)
    plot3(t1(1), t2(1), t3(1), 'o', ...
          'MarkerSize', 8, 'MarkerFaceColor', colors(t,:), ...
          'MarkerEdgeColor', colors(t,:), 'HandleVisibility','off');
    hold on

    % marker fine (triangolo)
    plot3(t1(end), t2(end), t3(end), '^', ...
          'MarkerSize', 8, 'MarkerFaceColor', colors(t,:), ...
          'MarkerEdgeColor', colors(t,:), 'HandleVisibility','off');

    xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
    ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
    zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
    grid on; axis equal;

end

h1 = plot3(nan, nan, nan, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
h2 = plot3(nan, nan, nan, '^', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% legend([h1 h2], {'Start', 'End'}, 'Location', 'northeastoutside');
