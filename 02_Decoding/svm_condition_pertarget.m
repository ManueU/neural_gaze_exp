clearvars -except mean_baseline std_baseline
% close all
clc


n_sets = 6;
n_arrays = 2;
n_trials = 32;
bin_size = 0.02;
period_pre = 0.1;
period_post = 0.5;
n_targets = 4; 

filename = {'../00_Data_extraction/free-gaze_BCI02_withtracker_exclUpdated.mat',... 
           '../00_Data_extraction/motor_BCI02_withtracker_exclUpdated.mat',...
           '../00_Data_extraction/controlled_BCI02_withtracker_exclUpdated.mat'};


%% Decoding condizione fissando il target
% Output: per ogni target -> acc, balAcc, cm

results = struct();
acc_t   = nan(n_targets,1);
bal_t   = nan(n_targets,1);
cm_t    = cell(n_targets,1);

labelsCond = ["Free-gaze","Gaze-on-center","Gaze-on-target"]; 

for tFix = 1:n_targets
    fprintf('\n=============================\n');
    fprintf('Condition decoding | fixed Target = %d\n', tFix);
    fprintf('=============================\n');

    X_all = [];
    Y_all = strings(0,1);

    for d = 1:numel(filename)
        fprintf('  Dataset: %s\n', filename{d});
        load(filename{d});

        [~, baseName, ext] = fileparts(filename{d});
        ds_name = [baseName ext];

        if strcmp(ds_name, 'controlled_BCI02_withtracker_exclUpdated.mat')
            PRE  = "Gaze";   POST = "Reach";
            cond_label = "Gaze-on-target";
        elseif strcmp(ds_name, 'motor_BCI02_withtracker_exclUpdated.mat')
            PRE  = "Pres12"; POST = "Reach";
            cond_label = "Gaze-on-center";
        else
            PRE  = "Pres12"; POST = "Reach";
            cond_label = "Free-gaze";
        end

        idx_pres  = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == PRE);
        idx_reach = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == POST);
        if isempty(idx_pres) || isempty(idx_reach)
            error('Stati PRE/POST non trovati: controlla PRE/POST');
        end

        n_bins_pre  = round(period_pre/bin_size);
        n_bins_post = round(period_post/bin_size);
        start_pre = size(data(1).Data(1).Interp(1).Task_states{idx_pres,2}, 1) - n_bins_pre + 1;
        end_post  = n_bins_post;

        % raccogli trial solo del target fissato
        for set = 1:n_sets
            for trial = 1:n_trials

                % usa array 1 per leggere Target_ID/Excluded (di solito identici per array)
                if data(set).Data(1).Interp(trial).Excluded ~= 0
                    continue
                end
                if data(set).Data(1).Interp(trial).Target_ID ~= tFix
                    continue
                end

                tmp_pre  = [];
                tmp_post = [];
                for array = 1:n_arrays
                    tmp_pre  = [tmp_pre,  data(set).Data(array).Interp(trial).Task_states{idx_pres,2}(start_pre:end, :)];
                    tmp_post = [tmp_post, data(set).Data(array).Interp(trial).Task_states{idx_reach,2}(1:end_post, :)];
                end

                matrix = [tmp_pre; tmp_post];       % [time x features]
                feat   = mean(matrix./bin_size, 1); % 1 x features

                X_all = [X_all; feat];
                Y_all = [Y_all; cond_label];
            end
        end
    end

    % controlla che ci siano campioni per tutte le condizioni
    Ycat = categorical(Y_all, labelsCond); % fissa ordine classi
    counts = countcats(Ycat);

    disp(table(categories(Ycat), counts, 'VariableNames', {'Condition','Count'}));

    if any(counts < 2)
        warning('Target %d: una o più condizioni hanno <2 campioni, salto.', tFix);
        continue
    end

    % === Bilanciamento per condizione (downsample) ===
    minCount = min(counts);

    rng(0);
    keep_idx = false(size(Ycat));
    for ci = 1:numel(labelsCond)
        idxC = find(Ycat == labelsCond(ci));
        idxC = idxC(randperm(numel(idxC), minCount));
        keep_idx(idxC) = true;
    end

    Xb = X_all(keep_idx,:);
    Yb = Ycat(keep_idx);

    % k-fold adattivo (non > minCount)
    k_fold = min(5, minCount);
    if k_fold < 2
        warning('Target %d: k_fold < 2 dopo bilanciamento, salto.', tFix);
        continue
    end

    fprintf('Balanced per condition: %d per classe | k_fold=%d\n', minCount, k_fold);

    % === SVM CV ===
    [acc, cm, metrics] = svm_cv(Xb, Yb, k_fold);

    acc_t(tFix) = acc;
    bal_t(tFix) = metrics.balancedAccuracy;
    cm_t{tFix}  = cm;

    fprintf('Target %d | Acc=%.1f%% | BalAcc=%.1f%% | MacroF1=%.2f\n', ...
        tFix, acc*100, metrics.balancedAccuracy*100, metrics.macroF1);

    % Plot confusion (opzionale)
    figure('Color','w');
    cc = confusionchart(cm, categories(Yb), 'Normalization','row-normalized', ...
        'RowSummary','off','ColumnSummary','off');
    cc.GridVisible = 'off';
    cc.XLabel = '\bf Predicted Condition';
    cc.YLabel = '\bf True Condition';
    title(sprintf('Fixed Target %d | BalAcc %.1f%%', tFix, metrics.balancedAccuracy*100));
end

results.acc = acc_t;
results.balAcc = bal_t;
results.cm = cm_t;

% Summary plot
figure('Color','w'); 
bar(100*bal_t); ylim([0 100]); grid on
xlabel('Target'); ylabel('Balanced Accuracy (%)');
title('Condition decoding (fixed target)');
yline(100/numel(labelsCond),'--','Chance');
