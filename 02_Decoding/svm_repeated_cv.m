function [results, cm_pooled] = svm_repeated_cv(X, Y, k_fold, n_repeats, seed)

    if nargin < 5
        seed = 0;
    end

    Y = categorical(Y);
    allClasses = categories(Y);
    nClasses = numel(allClasses);

    % Template SVM lineare
    t = templateSVM( ...
        'KernelFunction', 'linear', ...
        'Standardize', true);

    % Output per repeat
    acc_all      = nan(n_repeats,1);
    balAcc_all   = nan(n_repeats,1);
    macroF1_all  = nan(n_repeats,1);
    recall_all   = nan(n_repeats, nClasses);
    prec_all     = nan(n_repeats, nClasses);
    f1_all       = nan(n_repeats, nClasses);

    % Confusion matrix pooled su tutti i repeat
    cm_pooled = zeros(nClasses, nClasses);

    for r = 1:n_repeats

        rng(seed + r);
        cv = cvpartition(Y, 'KFold', k_fold);

        % categorici vuoti ma con le stesse categorie di Y
        Y_true_rep = Y([]);
        Y_pred_rep = Y([]);

        for k = 1:k_fold
            trainIdx = training(cv, k);
            testIdx  = test(cv, k);

            Ytr = Y(trainIdx);

            % Controllo classi presenti nel training
            missing = setdiff(allClasses, categories(removecats(Ytr)));
            if ~isempty(missing)
                error('Repeat %d, fold %d: classi mancanti nel TRAIN: %s', ...
                    r, k, strjoin(cellstr(missing), ', '));
            end

            % Pesi inversamente proporzionali alla frequenza
            classes_tr = categories(Ytr);
            counts_tr  = countcats(Ytr);
            Ntr        = numel(Ytr);
            Ktr        = numel(classes_tr);

            classW = Ntr ./ (Ktr * counts_tr);

            w = zeros(sum(trainIdx),1);
            for i = 1:Ktr
                w(Ytr == classes_tr{i}) = classW(i);
            end

            model = fitcecoc( ...
                X(trainIdx,:), Y(trainIdx), ...
                'Learners', t, ...
                'Coding', 'onevsall', ...
                'Weights', w);

            yhat = predict(model, X(testIdx,:));

            Y_true_rep = [Y_true_rep; Y(testIdx)];
            Y_pred_rep = [Y_pred_rep; yhat];
        end

        % Confusion matrix del repeat
        cm_rep = confusionmat(Y_true_rep, Y_pred_rep, 'Order', allClasses);
        cm_pooled = cm_pooled + cm_rep;

        % Metriche
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