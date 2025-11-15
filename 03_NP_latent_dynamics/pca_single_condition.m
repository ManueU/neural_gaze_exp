% =========================================================
% DESCRIZIONE:
% Esegue PCA sulle singole condizioni per verificare la separabilità
% nello spazio delle PC in base al target del movimento. 
%
% La finestra temporale analizzata è fissa rispetto all’onset del reach:
% w = [-100, +500] ms.
%
% Prima dell’esecuzione, modificare i parametri:
%   filename : nome del file .mat da caricare
%   PRE, POST: etichette delle condizioni di interesse
% =========================================================

clearvars 
close all
clc

sets_pca = [1,2,4,5,6];
n_sets = 5; 
n_arrays = 2;      
n_channels = 96;    
n_targets = 8; 
bin_size = 0.02; 
period_pre = 0.1; 
period_post = 0.5; 
PRE = "Pres12";
POST = "Reach";

%% Load 
filename = 'free-gaze_BCI02.mat';
load(filename)

%% Costruzione matrice PCA 
idx_pres = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == PRE); 
idx_reach = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == POST); 
if isempty(idx_pres) || isempty(idx_reach)
    error('Stati PRE/POST non trovati: controlla PRE/POST');
end

pre_bins  = max(1, round(period_pre/bin_size));
post_bins = max(1, round(period_post/bin_size));

pca_matrix_array = []; 
for array = 1:n_arrays
    pca_matrix = []; 

    for channel = 1:n_channels
        firing_rate = []; 
        
        for target = 1:n_targets
            M_spikes = [];

            for set_all = 1:n_sets
                set = sets_pca(set_all);
                idx = find([data(set).Data(array).Resampled.Target_ID] == target); 

                for j = 1:length(idx)
                    tmp_pre = data(set).Data(array).Resampled(idx(j)).Task_states{idx_pres, 2}(end-pre_bins+1:end, channel); 
                    tmp_post = data(set).Data(array).Resampled(idx(j)).Task_states{idx_reach,2}(1:post_bins, channel); 
                    matrix = [tmp_pre; tmp_post];  
                    
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


%% PCA
[Xz, muZ, sigmaZ] = zscore(pca_matrix_array, 0, 1); 
[coeff, score, latent, tsquared, explained] = pca(Xz, 'Algorithm','svd');


%% Figure (1) - PCA
nbin = pre_bins + post_bins;           

figure('Color','White'); hold on
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

w = 11; 
for t = 1:n_targets
    idx = (t-1)*nbin + (1:nbin);   
    traj = score(idx,1:3);        
    t1 = smoothdata(traj(:,1), 'gaussian', w);
    t2 = smoothdata(traj(:,2), 'gaussian', w);
    t3 = smoothdata(traj(:,3), 'gaussian', w);

    plot3(t1, t2, t3, ...
          'Color', colors(t,:), 'LineWidth', 2, 'MarkerFaceColor', colors(t,:));

    % marker inizio (cerchio)
    plot3(t1(1), t2(1), t3(1), 'o', ...
          'MarkerSize', 8, 'MarkerFaceColor', colors(t,:), ...
          'MarkerEdgeColor', colors(t,:), 'HandleVisibility','off');

    % marker fine (triangolo)
    plot3(t1(end), t2(end), t3(end), '^', ...
          'MarkerSize', 8, 'MarkerFaceColor', colors(t,:), ...
          'MarkerEdgeColor', colors(t,:), 'HandleVisibility','off');
end

h1 = plot3(nan, nan, nan, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
h2 = plot3(nan, nan, nan, '^', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
legend([h1 h2], {'Start', 'End'}, 'Location', 'northeastoutside');
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
grid on; axis equal;


%% Figure (2) - Loadings
nPC = 3;
figure('Color','White');
bar(abs(coeff(:,1:nPC)), 'stacked');   
xlabel('Channel'); ylabel('|loading|');
legend(compose('PC%d',1:nPC), 'Location','northeastoutside');
title('Loadings');
grid on;


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