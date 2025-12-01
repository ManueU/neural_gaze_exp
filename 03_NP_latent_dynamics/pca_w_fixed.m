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

clearvars
close all
clc

sets = [2,4,5,6]; 
n_sets = numel(sets); 
n_arrays = 2;      
n_channels = 96;    
n_targets = 8; 
bin_size = 0.02; 
period_pre = 1.0; 
period_post = 0.2; 

filename = {'../00_Data_extraction/free-gaze_BCI02.mat', ...
            '../00_Data_extraction/motor_BCI02.mat',...
            '../00_Data_extraction/controlled_BCI02.mat'};


%% Costruzione matrice PCA 
condition = []; 
for d = 1:numel(filename) 
    fprintf('\nDataset: %s\n', filename{d}); 
    load(filename{d});

    [~, baseName, ext] = fileparts(filename{d});
    ds_name = [baseName ext];

    if strcmp(ds_name, 'controlled_BCI02.mat')
        PRE = "Gaze";
        POST = "Reach";
    elseif strcmp(ds_name, 'gaze_BCI02.mat')
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
            firing_rate = []; 
            
            for target = 1:n_targets
                M_spikes = [];
    
                for set_ = 1:n_sets
                    set = sets(set_);
                    idx = find([data(set).Data(array).Interp.Target_ID] == target); 
    
                    for j = 1:length(idx)
                        tmp_pre = data(set).Data(array).Interp(idx(j)).Task_states{idx_pres, 2}(end-pre_bins+1:end, channel); 
                        tmp_post = data(set).Data(array).Interp(idx(j)).Task_states{idx_reach,2}(1:post_bins, channel); 
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
    condition = [condition; pca_matrix_array];
    condition_sep(d) = {pca_matrix_array};
end 

%% Filtraggio neuroni poco informativi
% varianza e firing rate medio per neurone
meanFR = mean(condition, 1);          % Hz
varFR  = var(condition, 0, 1);

minMeanFR = 1;       % Hz: neuroni con FR medio < 1 Hz vengono scartati
minVarFR  = 0.5;     % varianza minima (su firing rate in Hz)

valid_cols = (meanFR > minMeanFR) & (varFR > minVarFR);
fprintf('\nNeuroni totali: %d, neuroni tenuti dopo filtro: %d\n', numel(meanFR), sum(valid_cols));

condition = condition(:, valid_cols);
for d = 1:numel(condition_sep)
    condition_sep{d} = condition_sep{d}(:, valid_cols);
end

%% PCA
[Xz, muZ, sigmaZ] = zscore(condition, 0, 1); 
[coeff, score, latent, tsquared, explained] = pca(Xz, 'Algorithm','svd');

fprintf('\nVarianza spiegata:\n');
for k = [1 3 5 10 20 50]
    if k <= numel(explained)
        fprintf('  Prime %2d PC: %5.2f %%\n', k, sum(explained(1:k)));
    end
end
fprintf('  Totale     : %5.2f %%\n\n', sum(explained));


%% Figure (1) - PCA
nbin = pre_bins + post_bins;           
nCond = numel(filename);
T = repmat(n_targets * nbin, nCond, 1); 
offset = [0; cumsum(T(1:end-1))];

figure('Color','White'); hold on
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

w = 11; 
legend_entries = cell(1,4);
for c = 1:nCond
    for t = 1:n_targets
        rows = offset(c) + (t-1)*nbin + (1:nbin);
        traj = score(rows,1:3);
        t1 = smoothdata(traj(:,1), 'gaussian', w);
        t2 = smoothdata(traj(:,2), 'gaussian', w);
        t3 = smoothdata(traj(:,3), 'gaussian', w);

        plot3(t1, t2, t3, 'Color', colors(nCond,:), 'LineWidth', 2);
    
        % marker inizio (cerchio)
        plot3(t1(1), t2(1), t3(1), 'o', ...
              'MarkerSize', 8, 'MarkerFaceColor', colors(nCond,:), ...
              'MarkerEdgeColor', colors(t,:), 'HandleVisibility','off');
    
        % marker fine (triangolo)
        plot3(t1(end), t2(end), t3(end), '^', ...
              'MarkerSize', 8, 'MarkerFaceColor', colors(nCond,:), ...
              'MarkerEdgeColor', colors(t,:), 'HandleVisibility','off');
    end
end 

hStart = plot3(nan,nan,nan, 'o', 'Color','k', 'MarkerFaceColor','k');
hEnd   = plot3(nan,nan,nan, '^', 'Color','k', 'MarkerFaceColor','k');

if nCond == 1
    % Solo Start/End
    legend([hStart hEnd], {'Start','End'}, 'Location','northeastoutside');
% else
%     % Crea etichette leggibili dai filename
%     legEntries = cell(1, nCond + 2);  % condizioni + Start/End
%     styles = {'-','--',':','-.'};
%     hStyles = gobjects(1, nCond);
% 
%     for c = 1:nCond
%         % Estrai solo la parte prima del primo "_"
%         [~, baseName, ~] = fileparts(filename{c});
%         [namePart, ~] = strtok(baseName, '_');
%         legEntries{c} = namePart;
%         % Dummy per linea condizione
%         ls = styles{mod(c-1, numel(styles))+1};
%         hStyles(c) = plot3(nan,nan,nan, ls, 'Color','k', 'LineWidth', 2);
%     end
% 
%     legEntries{nCond+1} = 'Start';
%     legEntries{nCond+2} = 'End';
%     legend([hStyles hStart hEnd], legEntries, 'Location','northeastoutside');
end

xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
grid on; axis equal;


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






