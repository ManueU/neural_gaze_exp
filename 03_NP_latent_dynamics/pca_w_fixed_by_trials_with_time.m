% 1: medial arm
% 2: lateral hand

clearvars -except mean_baseline_common std_baseline_common
close all
clc

sets = [1,2,3,4,5,6];
n_sets = numel(sets);
n_arrays = 2;
n_targets = 8;
n_trials = 16;
bin_size = 0.02;
period_pre = 0.1;
period_post = 0.5;

filename = { ...
    '../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat', ...
    '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat', ...
    '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat'
    };

%% Definizione condizioni, marker e stati PRE/POST
nCond = numel(filename);

cond_labels = cell(nCond,1);
cond_markers = cell(nCond,1);
PRE_list = strings(nCond,1);
POST_list = strings(nCond,1);

for d = 1:nCond
    [~, baseName, ext] = fileparts(filename{d});
    ds_name = [baseName ext];
    ds_name_low = lower(ds_name);

    if contains(ds_name_low, 'free-gaze')
        cond_labels{d}  = 'Free-gaze';
        cond_markers{d} = 'o';   % tondo
        PRE_list(d)  = "Pres12";
        POST_list(d) = "Reach";

    elseif contains(ds_name_low, 'motor')
        cond_labels{d}  = 'Gaze-on-center';
        cond_markers{d} = 's';   % quadrato
        PRE_list(d)  = "Pres12";
        POST_list(d) = "Reach";

    elseif contains(ds_name_low, 'controlled')
        cond_labels{d}  = 'Gaze-on-target';
        cond_markers{d} = '^';   % triangolo
        PRE_list(d)  = "Gaze";
        POST_list(d) = "Reach";

    else
        error('Condizione non riconosciuta per il file: %s', ds_name);
    end
end

%% Costruzione matrice PCA
condition = [];
Y_all = [];
trial_cond_idx = [];   % salva a quale condizione appartiene ogni trial valido
valid_trials = zeros(nCond,1);

for d = 1:nCond
    pca_matrix = [];

    fprintf('\nDataset: %s\n', filename{d});
    load(filename{d});

    PRE  = PRE_list(d);
    POST = POST_list(d);

    Y = [];
    for set = 1:n_sets
        idx_valid = [data(set).Data(1).Interp.Excluded] == 0;
        Y = [Y; [data(set).Data(1).Interp(idx_valid).Target_ID]'];
    end

    valid_trials(d) = length(Y);
    Y_all = [Y_all; Y];
    trial_cond_idx = [trial_cond_idx; d*ones(length(Y),1)];

    idx_pres  = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == PRE);
    idx_reach = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == POST);

    if isempty(idx_pres) || isempty(idx_reach)
        error('Stati PRE/POST non trovati nel file %s', filename{d});
    end

    n_bins_pre  = round(period_pre/bin_size);
    n_bins_post = round(period_post/bin_size);

    start_pre = size(data(1).Data(1).Interp(1).Task_states{idx_pres,2}, 1) - n_bins_pre + 1;
    end_post  = n_bins_post;

    for set = 1:n_sets
        for trial = 1:n_trials
            idx_abs = (set-1)*n_trials + trial;
            disp(idx_abs)

            tmp_pre = [];
            tmp_post = [];
            matrix = [];

            if data(set).Data(1).Interp(trial).Excluded == 0
                for array = 1:n_arrays
                    tmp_pre  = [tmp_pre,  data(set).Data(array).Interp(trial).Task_states{idx_pres,2}(start_pre:end, :)];
                    tmp_post = [tmp_post, data(set).Data(array).Interp(trial).Task_states{idx_reach,2}(1:end_post, :)];
                end

                matrix = [tmp_pre; tmp_post];

                X = (mean(matrix./bin_size, 1) - mean_baseline_common)./std_baseline_common;
                X(isnan(X) | isinf(X)) = 0;

                pca_matrix = [pca_matrix; X];
            end
        end
    end

    condition = [condition; pca_matrix];
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

subplot(2,1,2)
plot(cumsum(explained), '-o','LineWidth',1.5)
xlabel('Number of Components')
ylabel('Cumulative Explained Variance (%)')
title('Cumulative Explained Variance')
xlim([1 size(score,2)])
ylim([0 100])
grid on

%% Preparazione per figure
% nBins = n_bins_pre + n_bins_post;
nBins = 1;
nTrials = sum(valid_trials);
nPC = 3;
S = score(:,1:nPC);

PC1 = reshape(S(:,1), [nBins, nTrials]);
PC2 = reshape(S(:,2), [nBins, nTrials]);
PC3 = reshape(S(:,3), [nBins, nTrials]);

w = 10;
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

%% Plot 1: traiettorie 3D di tutti i trial
figure('Color','w'); hold on; grid on;

for tr = 1:nTrials
    t1 = smoothdata(PC1(:,tr), 'gaussian', w);
    t2 = smoothdata(PC2(:,tr), 'gaussian', w);
    t3 = smoothdata(PC3(:,tr), 'gaussian', w);

    plot3(t1, t2, t3, '-', ...
        'LineWidth', 1, ...
        'Color', colors(Y_all(tr),:));

    % marker inizio
    plot3(t1(1), t2(1), t3(1), 'o', ...
        'MarkerSize', 7, ...
        'MarkerFaceColor', colors(Y_all(tr),:), ...
        'MarkerEdgeColor', colors(Y_all(tr),:), ...
        'HandleVisibility','off');

    % marker fine
    plot3(t1(end), t2(end), t3(end), '^', ...
        'MarkerSize', 7, ...
        'MarkerFaceColor', colors(Y_all(tr),:), ...
        'MarkerEdgeColor', colors(Y_all(tr),:), ...
        'HandleVisibility','off');
end

xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
% title('Traiettorie 3D dei trial')
view(45,25);

%% Plot 2: un pallino per trial con marker dipendente dalla condizione
figure('Color','w'); hold on; grid on;

for tr = 1:nTrials
    t1 = smoothdata(PC1(:,tr), 'gaussian', w);
    t2 = smoothdata(PC2(:,tr), 'gaussian', w);
    t3 = smoothdata(PC3(:,tr), 'gaussian', w);

    % centro medio della traiettoria
    p1 = mean(t1);
    p2 = mean(t2);
    p3 = mean(t3);

    this_cond = trial_cond_idx(tr);
    this_marker = cond_markers{this_cond};

    plot3(p1, p2, p3, this_marker, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', colors(Y_all(tr),:), ...
        'MarkerEdgeColor', colors(Y_all(tr),:));
end

xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
view(45,25);

% legenda condizioni
h_leg = gobjects(nCond,1);
for d = 1:nCond
    h_leg(d) = plot3(nan, nan, nan, cond_markers{d}, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'k');
end
legend(h_leg, cond_labels, 'Location', 'best');