clearvars -except mean_baseline_common std_baseline_common
close all
clc

n_sets = 6;
n_arrays = 2;
n_trials = 16;
bin_size = 0.02;
period_pre = 0.1;
period_post = 0.5;   % se 0, non usare feature POST

filename = { ...
    '../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat', ...
    '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat', ...
    '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat'
    % , ...
    % '../00_Data_extraction/BCI02_Session_0924/gaze_BCI02_exclUpdated.mat'
    };

%% Decoding condizione aggregando i trial di tutti i dataset
X_all = [];
Y_all = strings(0,1);

n_bins_pre  = round(period_pre / bin_size);
n_bins_post = round(period_post / bin_size);

for d = 1:numel(filename)
    fprintf('\nDataset: %s\n', filename{d});
    S = load(filename{d});
    data = S.data;

    [~, baseName, ext] = fileparts(filename{d});
    ds_name = [baseName ext];

    if strcmp(ds_name, 'controlled_BCI02_exclUpdated.mat')
        PRE  = "Gaze";
        POST = "Reach";
        cond_label = "Gaze-on-target";
    elseif strcmp(ds_name, 'motor_BCI02_exclUpdated.mat')
        PRE  = "Pres12";
        POST = "Reach";
        cond_label = "Gaze-on-center";
    elseif strcmp(ds_name, 'gaze_BCI02_exclUpdated.mat')
        PRE  = "Pres12";
        POST = "Gaze";
        cond_label = "Gaze-only";
    else
        PRE  = "Pres12";
        POST = "Reach";
        cond_label = "Free-gaze";
    end

    idx_pres  = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == PRE, 1);
    idx_reach = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == POST, 1);

    if isempty(idx_pres) || isempty(idx_reach)
        error('Stati PRE/POST non trovati nel dataset %s.', ds_name);
    end

    X_d = [];
    Y_d = strings(0,1);

    for set = 1:n_sets
        for trial = 1:n_trials

            if data(set).Data(1).Interp(trial).Excluded ~= 0
                continue;
            end

            tmp_pre = [];
            tmp_post = [];

            valid_trial = true;

            for array = 1:n_arrays
                pre_mat = data(set).Data(array).Interp(trial).Task_states{idx_pres,2};

                if size(pre_mat,1) < n_bins_pre
                    warning('Dataset %s | set %d trial %d array %d: PRE troppo corto, trial saltato.', ...
                        ds_name, set, trial, array);
                    valid_trial = false;
                    break;
                end

                start_pre = size(pre_mat,1) - n_bins_pre + 1;
                tmp_pre = [tmp_pre, pre_mat(start_pre:end, :)];

                if n_bins_post > 0
                    post_mat = data(set).Data(array).Interp(trial).Task_states{idx_reach,2};

                    if size(post_mat,1) < n_bins_post
                        warning('Dataset %s | set %d trial %d array %d: POST troppo corto, trial saltato.', ...
                            ds_name, set, trial, array);
                        valid_trial = false;
                        break;
                    end

                    tmp_post = [tmp_post, post_mat(1:n_bins_post, :)];
                end
            end

            if ~valid_trial
                continue;
            end

            mu = mean_baseline_common(:)';
            sigma = std_baseline_common(:)';
            sigma(sigma == 0) = 1;
            
            if n_bins_post > 0
                tmp_all = [tmp_pre; tmp_post];
            
                % Converto in firing rate e faccio una sola media temporale
                feat = mean(tmp_all ./ bin_size, 1);
                feat = (feat - mu) ./ sigma;
            else
                % Solo PRE
                feat = mean(tmp_pre ./ bin_size, 1);            
                feat = (feat - mu) ./ sigma;
            end

            if any(isnan(feat)) || any(isinf(feat))
                warning('Dataset %s | set %d trial %d: feature non valide, trial saltato.', ...
                    ds_name, set, trial);
                continue;
            end

            X_d = [X_d; feat];
            Y_d = [Y_d; cond_label];
        end
    end

    X_all = [X_all; X_d];
    Y_all = [Y_all; Y_d];
end

% Controllo finale
if isempty(X_all) || isempty(Y_all)
    error('Nessun campione valido trovato.');
end

Ycat = categorical(Y_all);
counts = countcats(Ycat);
minCount = min(counts);

k_fold = min(4, minCount);
if k_fold < 2
    error('Una o più condizioni hanno <2 campioni: impossibile fare k-fold CV.');
end

fprintf('\nCondition counts:\n');
disp(table(categories(Ycat), counts, 'VariableNames', {'Condition','Count'}));
fprintf('Using k_fold = %d (min class count = %d)\n', k_fold, minCount);

rng(0);
[acc, cm, metrics] = svm_cv_percond(X_all, Ycat, k_fold);

figure('Color','w');
classes = categories(Ycat);

cc = confusionchart(cm, classes, ...
    'Normalization', 'row-normalized', ...
    'RowSummary', 'off', ...
    'ColumnSummary', 'off');

cc.GridVisible = 'off';
cc.DiagonalColor = [0 0.35 0.7];
cc.OffDiagonalColor = [0.8 0.8 0.8];
cc.XLabel = '\bf Predicted Condition';
cc.YLabel = '\bf True Condition';

title(sprintf('Cond decoding | Acc: %.1f%% | BalAcc: %.1f%%', ...
    acc*100, metrics.balancedAccuracy*100));