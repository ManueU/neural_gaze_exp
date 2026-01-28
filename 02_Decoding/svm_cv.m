function [acc, cm, metrics] = svm_cv(X, Y, k_fold)
    Y = categorical(Y);
    cv = cvpartition(Y, 'KFold', k_fold);
    Y_true = categorical([]);
    Y_pred = categorical([]);

    t = templateSVM('KernelFunction','rbf','KernelScale','auto','Standardize',true);

    for k = 1:k_fold
        trainIdx = training(cv, k);
        testIdx  = test(cv, k);

        Ytr = Y(trainIdx);

        classes = categories(Ytr);
        counts  = countcats(Ytr);
        N       = numel(Ytr);
        K       = numel(classes);

        classW = N ./ (K * counts);

        w = zeros(size(Ytr));
        for i = 1:K
            w(Ytr == classes{i}) = classW(i);
        end

        model = fitcecoc(X(trainIdx,:), Y(trainIdx), 'Learners', t, 'Coding', 'onevsall', 'Weights', w);
        yhat  = predict(model, X(testIdx,:));

        Y_true = [Y_true; Y(testIdx)];
        Y_pred = [Y_pred; yhat];
    end

    acc = mean(Y_true == Y_pred);
    allClasses = categories(Y);
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

