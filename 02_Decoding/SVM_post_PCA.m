clearvars -except mean_baseline_common std_baseline_common 
clc
% close all

n_sets = 6;
n_arrays = 2;
n_trials = 16;
bin_size = 0.02;
period_pre = 0.1;
period_post = 0.5;

filename = {'../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat',... 
            '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat',...
            '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat',...
            '../00_Data_extraction/BCI02_Session_0924/gaze_BCI02_exclUpdated.mat'};

acc_raw_all = cell(numel(filename),1);
acc_pca_all = cell(numel(filename),1);
metrics_raw_all = cell(numel(filename),1);
metrics_pca_all = cell(numel(filename),1);
cm_raw_all = cell(numel(filename),1);
cm_pca_all = cell(numel(filename),1);

n_repeats = 20;
seed = 0;
nPC = 20;   % da variare/testare

for d = 1:numel(filename)

    fprintf('\nDataset: %s\n', filename{d});
    load(filename{d});

    [~, baseName, ext] = fileparts(filename{d});
    ds_name = [baseName ext];

    if strcmp(ds_name, 'controlled_BCI02_exclUpdated.mat')
        PRE  = "Gaze";
        POST = "Reach";
    elseif strcmp(ds_name, 'gaze_BCI02_exclUpdated.mat')
        PRE  = "Pres12";
        POST = "Gaze";
    else
        PRE  = "Pres12";
        POST = "Reach";
    end

    % -------------------------
    % Costruzione Y
    % -------------------------
    Y = [];
    for set = 1:n_sets
        idx = [data(set).Data(1).Interp.Excluded] == 0;
        Y = [Y; [data(set).Data(1).Interp(idx).Target_ID]'];
    end

    % -------------------------
    % Indici stati
    % -------------------------
    idx_pre  = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == PRE);
    idx_post = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == POST);

    if isempty(idx_pre) || isempty(idx_post)
        error('Stati PRE/POST non trovati.');
    end

    n_bins_pre  = round(period_pre/bin_size);
    n_bins_post = round(period_post/bin_size);

    % numero trial validi
    n_valid = sum(arrayfun(@(s) sum([data(s).Data(1).Interp.Excluded] == 0), 1:n_sets));

    X = cell(n_valid,1);
    j = 1;

    for set = 1:n_sets
        for trial = 1:n_trials

            if data(set).Data(1).Interp(trial).Excluded == 0

                tmp_pre = [];
                tmp_post = [];

                for array = 1:n_arrays
                    mat_pre = data(set).Data(array).Interp(trial).Task_states{idx_pre,2};
                    mat_post = data(set).Data(array).Interp(trial).Task_states{idx_post,2};

                    start_pre = size(mat_pre,1) - n_bins_pre + 1;

                    tmp_pre  = [tmp_pre,  mat_pre(start_pre:end, :)];
                    tmp_post = [tmp_post, mat_post(1:n_bins_post, :)];
                end

                matrix = [tmp_pre; tmp_post];              % bins x 192
                X{j} = mean(matrix./bin_size, 1);          % 1 x 192 (Hz)
                j = j + 1;
            end
        end
    end

    X = cell2mat(X);   % nTrialsValid x 192
    X_norm = (X - mean_baseline_common) ./ std_baseline_common;
    X_norm(isnan(X_norm) | isinf(X_norm)) = 0;

    % -------------------------
    % k-fold adattivo
    % -------------------------
    Ycat = categorical(Y);
    counts = countcats(Ycat);
    minCount = min(counts);
    k_fold = min(6, minCount);

    if k_fold < 2
        error('Una o più classi hanno meno di 2 campioni.');
    end

    fprintf('Class counts: ');
    fprintf('%d ', counts);
    fprintf('\nUsing repeated %d-fold CV with %d repeats\n', k_fold, n_repeats);

    % -------------------------
    % RAW SVM
    % -------------------------
    [results_raw, cm_raw] = svm_repeated_cv(X_norm, Y, k_fold, n_repeats, seed);

    fprintf('RAW accuracy: %.3f ± %.3f\n', mean(results_raw.acc), std(results_raw.acc));
    fprintf('RAW balanced accuracy: %.3f ± %.3f\n', mean(results_raw.balAcc), std(results_raw.balAcc));

    % -------------------------
    % PCA + SVM
    % -------------------------
    [results_pca, cm_pca] = svm_pca_repeated_cv(X_norm, Y, k_fold, n_repeats, seed, nPC);

    fprintf('PCA accuracy: %.3f ± %.3f\n', mean(results_pca.acc), std(results_pca.acc));
    fprintf('PCA balanced accuracy: %.3f ± %.3f\n', mean(results_pca.balAcc), std(results_pca.balAcc));

    acc_raw_all{d} = results_raw.acc;
    acc_pca_all{d} = results_pca.acc;
    metrics_raw_all{d} = results_raw;
    metrics_pca_all{d} = results_pca;
    cm_raw_all{d} = cm_raw;
    cm_pca_all{d} = cm_pca;
end

%% Figure
figure('Color','w')
hold on

nCond = numel(metrics_raw_all);

% colori
col_raw = [0.3 0.5 0.8];
col_pca = [0.8 0.4 0.4];

bar_width = 0.35;

bal_raw_mean = zeros(nCond,1);
bal_raw_std  = zeros(nCond,1);

bal_pca_mean = zeros(nCond,1);
bal_pca_std  = zeros(nCond,1);

for d = 1:nCond
    
    % RAW
    bal_raw = metrics_raw_all{d}.balAcc * 100;
    bal_raw_mean(d) = mean(bal_raw);
    bal_raw_std(d)  = std(bal_raw);
    
    % PCA
    bal_pca = metrics_pca_all{d}.balAcc * 100;
    bal_pca_mean(d) = mean(bal_pca);
    bal_pca_std(d)  = std(bal_pca);
end

x = 1:nCond;

% barre
b1 = bar(x - bar_width/2, bal_raw_mean, bar_width, ...
    'FaceColor', col_raw, 'EdgeColor','none');

b2 = bar(x + bar_width/2, bal_pca_mean, bar_width, ...
    'FaceColor', col_pca, 'EdgeColor','none');

% errorbar
errorbar(x - bar_width/2, bal_raw_mean, bal_raw_std, ...
    'k','LineStyle','none','LineWidth',1.2)

errorbar(x + bar_width/2, bal_pca_mean, bal_pca_std, ...
    'k','LineStyle','none','LineWidth',1.2)

% scatter dei repeat
for d = 1:nCond
    
    % RAW
    y_raw = metrics_raw_all{d}.balAcc * 100;
    x_raw = (x(d) - bar_width/2) + (rand(size(y_raw))-0.5)*0.1;
    
    scatter(x_raw, y_raw, 20, ...
        'MarkerFaceColor', col_raw, ...
        'MarkerEdgeColor','none', ...
        'MarkerFaceAlpha',0.4)
    
    % PCA
    y_pca = metrics_pca_all{d}.balAcc * 100;
    x_pca = (x(d) + bar_width/2) + (rand(size(y_pca))-0.5)*0.1;
    
    scatter(x_pca, y_pca, 20, ...
        'MarkerFaceColor', col_pca, ...
        'MarkerEdgeColor','none', ...
        'MarkerFaceAlpha',0.4)
end

% chance level
nClassi = numel(unique(Y));
chance = (1/nClassi)*100;

yline(chance,'--k','Chance', ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','bottom')

ylim([0 100])
ylabel('Balanced accuracy (%)')
xlabel('Condition')

xticks(x)
xticklabels({'Free-gaze','Gaze-on-center','Gaze-on-target','Gaze-only'})

legend({'Raw (192 ch)','PCA + SVM'}, 'Location','northwest')

title('Decoding performance: Raw vs PCA')

grid on
box off