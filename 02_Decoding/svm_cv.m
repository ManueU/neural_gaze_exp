function [acc, cm, metrics] = svm_cv(X, Y, k_fold)
    Y = categorical(Y);
    allClasses = categories(Y);
    
    cv = cvpartition(Y, 'KFold', k_fold);
    Y_true = categorical([]);
    Y_pred = categorical([]);

    t = templateSVM('KernelFunction','linear','KernelScale','auto','Standardize',true);

    for k = 1:k_fold
        trainIdx = training(cv, k);
        testIdx  = test(cv, k);

        Ytr = Y(trainIdx);

        missing = setdiff(allClasses, categories(removecats(Ytr)));
        if ~isempty(missing)
            warning('Fold %d: classi mancanti nel TRAIN: %s', k, strjoin(cellstr(missing), ', '));
            error('Riduci k_fold o bilancia i dati: fold con classi mancanti.');
        end

        classes = categories(Ytr);
        counts  = countcats(Ytr);
        N       = numel(Ytr);
        K       = numel(classes);

        classW = N ./ (K * counts);

        w = zeros(sum(trainIdx),1);
        for i = 1:K
            w(Ytr == classes{i}) = classW(i);
        end

        model = fitcecoc(X(trainIdx,:), Y(trainIdx), 'Learners', t, 'Coding', 'onevsone', 'Weights', w);
        yhat  = predict(model, X(testIdx,:));

        Y_true = [Y_true; Y(testIdx)];
        Y_pred = [Y_pred; yhat];
    end

    acc = mean(Y_true == Y_pred);
    cm = confusionmat(Y_true, Y_pred, 'Order', allClasses);

    % Balanced accuracy (media dei recall per classe)
    tp = diag(cm);
    fn = sum(cm, 2) - tp;
    recall = tp ./ max(tp + fn, 1); % evita divisione per zero
    balAcc = mean(recall);

    % Macro-F1 (media F1 per classe)
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

