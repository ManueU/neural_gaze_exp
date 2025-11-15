function acc_dist = svm_resampled_cv(X, Y, k_fold, n_resamp)
% Esegue SVM one-vs-all con k-fold cross-validation ripetuta n_resamp volte
% Restituisce distribuzione di accuratezze

    acc_dist = zeros(n_resamp,1);
    for r = 1:n_resamp
        cv = cvpartition(Y, 'KFold', k_fold);
        acc_fold = zeros(k_fold,1);

        for k = 1:k_fold
            trainIdx = training(cv, k);
            testIdx  = test(cv, k);

            % Template SVM base
            t = templateSVM('KernelFunction','linear', 'Standardize',true);

            % SVM multiclasse one-vs-all
            model = fitcecoc(X(trainIdx,:), Y(trainIdx), ...
                             'Learners', t, ...
                             'Coding','onevsall');  

            % Predizione
            Ypred = predict(model, X(testIdx,:));

            acc_fold(k) = mean(Ypred == Y(testIdx));
        end

        acc_dist(r) = mean(acc_fold);
    end
end