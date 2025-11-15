function [acc, cm] = svm_cv(X, Y, k_fold)
    rng(99)    
    cv = cvpartition(Y, 'KFold', k_fold);
    Y_true = []; Y_pred = [];

    t = templateSVM('KernelFunction','rbf','KernelScale','auto','Standardize',true);
    for k = 1:k_fold
        trainIdx = training(cv, k);
        testIdx  = test(cv, k);

        model = fitcecoc(X(trainIdx,:), Y(trainIdx), 'Learners', t, 'Coding', 'onevsall');
        yhat  = predict(model, X(testIdx,:));

        Y_true = [Y_true; Y(testIdx)];
        Y_pred = [Y_pred; yhat];
    end

    acc = mean(Y_true == Y_pred);
    cm  = confusionmat(Y_true, Y_pred);
end