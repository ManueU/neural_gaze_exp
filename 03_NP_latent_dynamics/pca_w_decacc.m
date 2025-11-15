% =========================================================
% DESCRIZIONE:
% Esegue un’analisi PCA su più condizioni sperimentali per valutare:
% - la separabilità dei target di movimento nello spazio delle componenti principali (PC);
% - la relazione tra le diverse condizioni nello stesso spazio PC.
%
% La finestra temporale analizzata è calcolata sulla base del picco di
% accuratezza del decoder.
%
% Prima dell’esecuzione, modificare i parametri:
%   filename : nome del file .mat da caricare
%   PRE, POST: etichette delle condizioni di interesse
% =========================================================

clearvars
% close all
clc 

sets_decoder = [1, 2, 4, 5, 6];
sets_pca = [1, 2, 4, 5, 6];
n_sets = 5;
n_arrays = 2;  
n_channels = 96; 
n_trials = 32; 
n_targets = 8; 
bin_size = 0.02; 

filename = {'../00_Data_extraction/free-gaze_BCI02.mat', ...
            '../00_Data_extraction/motor_BCI02.mat'};
 
%% Decoding over time with SVM
for d = 1:numel(filename) 
    disp(filename(d)); 
    load(filename{d});

    % bins per trial
    N = size(data(1).Data(1).Resampled(1).Trial, 1); 
    rec_duration = N*bin_size; 

    % finestra scorrevole
    w_length = 0.6; 
    overlap = 0.5*w_length;
    N_w = round(w_length/bin_size);
    N_o = round(overlap/bin_size);
    
    % etichette Y
    Y = []; 
    for set_all = 1:n_sets
        set = sets_decoder(set_all);
        Y = [Y; [data(set).Data(1).Resampled.Target_ID]']; 
    end 
    classes = unique(Y,'stable');
    n_classes = numel(classes);
    
    % numero finestre
    n_acc = floor((N - N_w)/(N_w-N_o)) + 1; 
    acc_overall = zeros(n_acc,1);
    acc_class   = zeros(n_acc, n_classes);
    t = zeros(n_acc,1);
    cm = cell(n_acc,1);

    % loop sulle finestre
    start_w = 1; 
    end_w = start_w + N_w - 1; 
    
    for w = 1:n_acc
        X = cell(n_trials*n_sets,1);
        j = 1; 
        for set_all = 1:n_sets
            set = sets_decoder(set_all);
            for trial = 1:n_trials
                SVM_matrix = []; 
                for array = 1:n_arrays
                    SVM_matrix = [SVM_matrix, data(set).Data(array).Resampled(trial).Trial(start_w:end_w, :)]; 
                end 
                X{j} = mean(SVM_matrix,1)./bin_size;
                j = j + 1; 
            end   
        end 
        X = cell2mat(X); 
        
        % k-fold SVM 
        k_fold = 5; 
        [acc_overall(w), cm{w}] = svm_cv(X, Y, k_fold);

        % accuracy per classe (diag della cm normalizzata per riga) 
        cm_norm = cm{w} ./ max(sum(cm{w},2),1);
        acc_class(w,:) = diag(cm_norm)*100;
        
        t(w) = ((start_w + end_w - 1)/2) * bin_size;

        % prossima finestra
        start_w = start_w + (N_w - N_o); 
        end_w = start_w + N_w - 1; 
    end 

    % Figure
    figure('Color', 'White')
    events_time_tmp = []; 
    for i = 1:length(data(1).Data(2).Resampled(1).Task_states)
        events_time = [events_time_tmp; size(data(1).Data(2).Resampled(1).Task_states{i,2},1)*bin_size];
        events_time_tmp = events_time; 
    end 
    increment_times = cumsum(events_time); 
    labels = string(data(1).Data(2).Resampled(1).Task_states(:,1));

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
    alpha = 0.5;

    w_smooth = 5; 
    for c = 1:n_classes
        acc_smooth = smoothdata(acc_class(:, c), 'gaussian', w_smooth);
        plot(t', acc_smooth, 'LineWidth', 1.0, 'Color', [colors(c,:),  alpha], 'HandleVisibility','off'), hold on
    end
    acc_smooth_overall = smoothdata(acc_overall, 'gaussian', 4);
    plot(t', acc_overall*100, 'LineWidth', 1.5, 'Color', 'k', 'DisplayName', 'Overall'), hold on

    if exist('increment_times','var') && ~isempty(increment_times)
       xline(increment_times, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
       k = min(numel(labels), numel(increment_times));
       x_times = [0; increment_times(:)];
       x_text  = x_times(1:k) + diff(x_times(1:k+1))/2;
       ylim([0 100]);
       ax = gca; 
       y_text = (ax.YLim(2) - 10)*ones(1,length(x_text)); 
       text(x_text, y_text, labels(1:k), 'HorizontalAlignment','center');
    end

    yline((1/n_classes)*100,'-', 'Chance', 'HandleVisibility','off'); 
    legend show; 
    xlabel('Time (s)');
    ylabel('Accuracy (%)');
    xlim([0 rec_duration]);
    hold on

    % salva risultati
    acc_class_files{d} = acc_class;
    acc_overall_files{d} = [t, acc_overall]; 

    % picco per definire finestra (centro del massimo)
    [maxv, id_all] = max(acc_overall_files{d}(:,2));
    ids = find(acc_overall_files{d}(:,2) == maxv);
    
    if numel(ids) >= 2
        t_des = mean(acc_overall_files{d}(ids,1));
    else
        t_des = acc_overall_files{d}(id_all,1);
    end
    N_peak = round(t_des / bin_size);
    
    plot(t(ids), acc_overall_files{1,d}(ids,2)*100, 'r*', 'HandleVisibility','off'); 
    hold on

    delta_t = 0.25; % +- 250 ms
    t_pre(d)  = t_des - delta_t;                  
    t_post(d) = t_des + delta_t;
    xline(t_pre(d),  '--r', 'LineWidth', 1.0, ...
       'LabelOrientation','horizontal', 'HandleVisibility','off');
    xline(t_post(d), '--r', 'LineWidth', 1.0, ...
           'LabelOrientation','horizontal', 'HandleVisibility','off');
end

%% Costruzione matrice PCA 
condition = []; 
for d = 1:numel(filename) 
    fprintf('\nDataset: %s\n', filename{d}); 
    load(filename{d});

    start_idx = round(t_pre(d)/bin_size);
    end_idx = round(t_post(d)/bin_size);

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
                        matrix = data(set).Data(array).Resampled(idx(j)).Trial(start_idx:end_idx, channel); 
                        
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

col_std = std(condition, 0, 1);          % deviazione standard per colonna
valid_cols = col_std > 0;                % colonne con varianza > 0
condition = condition(:, valid_cols);
for d = 1:numel(condition_sep)
    condition_sep{d} = condition_sep{d}(:, valid_cols);
end

%% PCA
[Xz, muZ, sigmaZ] = zscore(condition, 0, 1); 
[coeff, score, latent, tsquared, explained] = pca(Xz, 'Algorithm','svd');


%% Figure (1) - PCA
nbin = end_idx - start_idx;           
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

        plot3(t1, t2, t3, 'Color', colors(t,:), 'LineWidth', 2, ...
              'LineStyle', styles{mod(c-1,numel(styles))+1});
    
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

hStart = plot3(nan,nan,nan, 'o', 'Color','k', 'MarkerFaceColor','k');
hEnd   = plot3(nan,nan,nan, '^', 'Color','k', 'MarkerFaceColor','k');

if nCond == 1
    % Solo Start/End
    legend([hStart hEnd], {'Start','End'}, 'Location','northeastoutside');
else
    % Crea etichette leggibili dai filename
    legEntries = cell(1, nCond + 2);  % condizioni + Start/End
    styles = {'-','--',':','-.'};
    hStyles = gobjects(1, nCond);

    for c = 1:nCond
        % Estrai solo la parte prima del primo "_"
        [namePart, ~] = strtok(filename{c}, '_');
        legEntries{c} = namePart;
        % Dummy per linea condizione
        hStyles(c) = plot3(nan,nan,nan, styles{c}, 'Color','k', 'LineWidth',2);
    end

    legEntries{nCond+1} = 'Start';
    legEntries{nCond+2} = 'End';
    legend([hStyles hStart hEnd], legEntries, 'Location','northeastoutside');
end

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


%% Saving video
% outputVideo = VideoWriter('motor_free_BCI02.mp4', 'Uncompressed AVI'); 
% outputVideo.FrameRate = 10;  
% open(outputVideo);
% 
% for angle = 0:2:360
%     view(angle, 30);   
%     frame = getframe(gcf);  
%     writeVideo(outputVideo, frame);  
% end
% close(outputVideo);

%% Procrustes alignment
% Z-score usando media e sigma globali (muZ, sigmaZ)
Xz_cond1 = (condition_sep{1} - muZ) ./ sigmaZ;
Xz_cond2 = (condition_sep{2} - muZ) ./ sigmaZ;

% Proiezione sulle PC comuni
scores1 = Xz_cond1 * coeff;   
scores2 = Xz_cond2 * coeff;

% Parametri per il confronto Procrustes
pc_idx = 1:3;                     % PC da usare 
T = nbin;   
K = numel(pc_idx);  

% Seleziona solo le PC desiderate
scores1 = scores1(:, pc_idx);     
scores2 = scores2(:, pc_idx);     

% Procrustes per ogni target
d_targets = zeros(n_targets,1);
for tg = 1:n_targets
    % righe di X/Y corrispondenti al target tg
    idx = (tg-1)*T + (1:T);    

    Xt = scores1(idx, :);     
    Yt = scores2(idx, :);      

    % Procrustes per il singolo target
    [d_t, Yt_aligned, transform_t] = procrustes(Xt, Yt, ...
        'scaling', true, 'reflection', false);

    d_targets(tg) = d_t;
end
disp('Distanza di Procrustes per target:');
disp(d_targets);

%% Figure (2)
figure('Color', 'w')
bar(d_targets);
xlabel('Target');
ylabel('Distanza di Procrustes');
title('Similarità della forma della traiettoria tra condizioni');
