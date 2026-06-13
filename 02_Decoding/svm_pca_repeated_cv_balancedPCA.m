function [results, cm_pooled] = svm_pca_repeated_cv_balancedPCA(X, Y, k_fold, n_repeats, seed, nPC)

    if nargin < 5
        seed = 0;
    end

    Y = categorical(Y(:));
    allClasses = categories(Y);
    nClasses = numel(allClasses);

    % Output per repeat
    acc_all      = nan(n_repeats,1);
    balAcc_all   = nan(n_repeats,1);
    macroF1_all  = nan(n_repeats,1);
    recall_all   = nan(n_repeats, nClasses);
    prec_all     = nan(n_repeats, nClasses);
    f1_all       = nan(n_repeats, nClasses);

    % Confusion matrix pooled
    cm_pooled = zeros(nClasses, nClasses);

    for r = 1:n_repeats

        rng(seed + r);
        cv = cvpartition(Y, 'KFold', k_fold, 'Stratify', true);

        Y_true_rep = Y([]);
        Y_pred_rep = Y([]);

        for f = 1:k_fold

            trainIdx = training(cv, f);
            testIdx  = test(cv, f);

            Xtr = X(trainIdx, :);
            Xte = X(testIdx, :);
            Ytr = Y(trainIdx);
            Yte = Y(testIdx);

            % Controllo classi presenti nel training
            missing = setdiff(allClasses, categories(removecats(Ytr)));
            if ~isempty(missing)
                error('Repeat %d, fold %d: classi mancanti nel TRAIN: %s', ...
                    r, f, strjoin(cellstr(missing), ', '));
            end

            % =========================================================
            % 1) PCA stimata su subset bilanciato del training set
            % =========================================================
            classes_tr = categories(Ytr);
            counts_tr  = countcats(Ytr);
            minCount   = min(counts_tr);

            idx_bal = [];

            for i = 1:numel(classes_tr)
                idx_i = find(Ytr == classes_tr{i});

                % campionamento casuale senza replacement
                idx_i = idx_i(randperm(numel(idx_i), minCount));

                idx_bal = [idx_bal; idx_i(:)];
            end

            % mescola gli indici bilanciati
            idx_bal = idx_bal(randperm(numel(idx_bal)));

            Xtr_bal = Xtr(idx_bal, :);
            Ytr_bal = Ytr(idx_bal); %#ok<NASGU>

            % PCA solo sul training bilanciato
            [coeff, ~, ~, ~, ~, mu_pca] = pca(Xtr_bal, 'Algorithm', 'svd');

            k = min(nPC, size(coeff,2));
            if k < 1
                error('Repeat %d, fold %d: numero di PC non valido.', r, f);
            end

            % =========================================================
            % 2) Proiezione di TUTTO il training e del test
            %    nello spazio PCA stimato sul training bilanciato
            % =========================================================
            Xtr_pca = (Xtr - mu_pca) * coeff(:,1:k);
            Xte_pca = (Xte - mu_pca) * coeff(:,1:k);

            % =========================================================
            % 3) Pesi per SVM su tutto il training set
            % =========================================================
            counts_tr_full = countcats(Ytr);
            Ntr_full       = numel(Ytr);
            Ktr_full       = numel(classes_tr);

            classW = Ntr_full ./ (Ktr_full * counts_tr_full);

            w = zeros(Ntr_full,1);
            for i = 1:Ktr_full
                w(Ytr == classes_tr{i}) = classW(i);
            end

            % =========================================================
            % 4) SVM multiclass
            % =========================================================
            t = templateSVM( ...
                'KernelFunction', 'linear', ...
                'Standardize', false);

            model = fitcecoc( ...
                Xtr_pca, Ytr, ...
                'Learners', t, ...
                'Coding', 'onevsall', ...
                'Weights', w);

            yhat = predict(model, Xte_pca);

            Y_true_rep = [Y_true_rep; Yte];
            Y_pred_rep = [Y_pred_rep; yhat];
        end

        % =============================================================
        % Metriche sul repeat
        % =============================================================
        cm_rep = confusionmat(Y_true_rep, Y_pred_rep, 'Order', allClasses);
        cm_pooled = cm_pooled + cm_rep;

        acc_rep = mean(Y_true_rep == Y_pred_rep);

        tp = diag(cm_rep);
        fn = sum(cm_rep, 2) - tp;
        fp = sum(cm_rep, 1)' - tp;

        recall = tp ./ max(tp + fn, 1);
        precision = tp ./ max(tp + fp, 1);
        f1 = 2 * (precision .* recall) ./ max(precision + recall, eps);

        balAcc = mean(recall);
        macroF1 = mean(f1);

        acc_all(r) = acc_rep;
        balAcc_all(r) = balAcc;
        macroF1_all(r) = macroF1;
        recall_all(r,:) = recall';
        prec_all(r,:) = precision';
        f1_all(r,:) = f1';
    end

    results = struct();
    results.acc = acc_all;
    results.balAcc = balAcc_all;
    results.macroF1 = macroF1_all;
    results.perClassRecall = recall_all;
    results.perClassPrecision = prec_all;
    results.perClassF1 = f1_all;
    results.classOrder = allClasses;
end