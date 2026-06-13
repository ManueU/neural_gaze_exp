function [acc, cm, metrics] = svm_cv_percond(X, Y, k_fold)

    Y = categorical(Y);
    allClasses = categories(Y);

    cv = cvpartition(Y, 'KFold', k_fold);

    Y_true = categorical(strings(0,1), allClasses);
    Y_pred = categorical(strings(0,1), allClasses);

    t = templateSVM( ...
        'KernelFunction', 'linear', ...
        'KernelScale', 'auto', ...
        'Standardize', false);

    for k = 1:k_fold
        trainIdx = training(cv, k);
        testIdx  = test(cv, k);

        Ytr = Y(trainIdx);
        presentClasses = categories(removecats(Ytr));
        missing = setdiff(allClasses, presentClasses);

        if ~isempty(missing)
            warning('Fold %d: classi mancanti nel TRAIN: %s', ...
                k, strjoin(cellstr(missing), ', '));
            error('Riduci k_fold o bilancia i dati: fold con classi mancanti.');
        end

        classes = categories(removecats(Ytr));
        counts = countcats(removecats(Ytr));
        N = numel(Ytr);
        K = numel(classes);

        classW = N ./ (K * counts);

        w = zeros(sum(trainIdx), 1);
        for i = 1:K
            w(Ytr == classes{i}) = classW(i);
        end

        model = fitcecoc( ...
            X(trainIdx,:), Y(trainIdx), ...
            'Learners', t, ...
            'Coding', 'onevsone', ...
            'Weights', w);

        yhat = predict(model, X(testIdx,:));

        Y_true = [Y_true; categorical(Y(testIdx), allClasses)];
        Y_pred = [Y_pred; categorical(yhat, allClasses)];
    end

    acc = mean(Y_true == Y_pred);
    cm = confusionmat(Y_true, Y_pred, 'Order', categorical(allClasses));

    % Balanced accuracy = media dei recall per classe
    tp = diag(cm);
    fn = sum(cm, 2) - tp;
    recall = tp ./ max(tp + fn, 1);
    balAcc = mean(recall);

    % Macro-F1
    fp = sum(cm, 1)' - tp;
    precision = tp ./ max(tp + fp, 1);
    f1 = 2 * (precision .* recall) ./ max(precision + recall, eps);
    macroF1 = mean(f1);

    metrics = struct( ...
        'balancedAccuracy', balAcc, ...
        'macroF1', macroF1, ...
        'perClassRecall', recall, ...
        'perClassPrecision', precision, ...
        'perClassF1', f1, ...
        'classOrder', {allClasses});
end