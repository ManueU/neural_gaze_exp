function [results, cm_pooled] = svm_pca_repeated_cv(X, Y, k_fold, n_repeats, seed, nPC)

    rng(seed);

    Y = Y(:);
    classes = unique(Y);
    nClass = numel(classes);

    acc = zeros(n_repeats,1);
    balAcc = zeros(n_repeats,1);
    macroF1 = zeros(n_repeats,1);

    cm_pooled = zeros(nClass, nClass);

    for r = 1:n_repeats

        rng(seed + r);
        cv = cvpartition(Y, 'KFold', k_fold, 'Stratify', true);

        y_true_all = [];
        y_pred_all = [];

        for f = 1:k_fold

            idxTrain = training(cv, f);
            idxTest  = test(cv, f);

            Xtr = X(idxTrain, :);
            Xte = X(idxTest, :);
            Ytr = Y(idxTrain);
            Yte = Y(idxTest);

            Ytr_cat = categorical(Ytr);
            classes_tr = categories(Ytr_cat);
            counts_tr  = countcats(Ytr_cat);
            Ntr        = numel(Ytr_cat);
            Ktr        = numel(classes_tr);
            
            classW = Ntr ./ (Ktr * counts_tr);
            
            w = zeros(Ntr,1);
            for i = 1:Ktr
                w(Ytr_cat == classes_tr{i}) = classW(i);
            end

            % ---------------------------------
            % PCA SOLO sul training set
            % ---------------------------------
            [coeff, score_tr, ~, ~, ~, mu_pca] = pca(Xtr, 'Algorithm', 'svd');
         
            k = min(nPC, size(score_tr,2));

            Xtr_pca = score_tr(:,1:k);
            Xte_pca = (Xte - mu_pca) * coeff(:,1:k);

            % ---------------------------------
            % SVM multiclass
            % ---------------------------------
            t = templateSVM('KernelFunction','linear', ...
                            'Standardize',false);

            mdl = fitcecoc(Xtr_pca, Ytr, ...
                           'Learners', t, ...
                           'Coding', 'onevsall', ...
                           'Weights', w);

            Ypred = predict(mdl, Xte_pca);

            y_true_all = [y_true_all; Yte];
            y_pred_all = [y_pred_all; Ypred];
        end

        % confusion matrix repeat
        cm = confusionmat(y_true_all, y_pred_all, 'Order', classes);
        cm_pooled = cm_pooled + cm;

        % metriche
        acc(r) = sum(diag(cm)) / sum(cm(:));

        recall = diag(cm) ./ sum(cm,2);
        recall(isnan(recall)) = 0;
        balAcc(r) = mean(recall);

        precision = diag(cm) ./ sum(cm,1)';
        precision(isnan(precision)) = 0;

        f1 = 2 * (precision .* recall) ./ (precision + recall);
        f1(isnan(f1)) = 0;
        macroF1(r) = mean(f1);
    end

    results.acc = acc;
    results.balAcc = balAcc;
    results.macroF1 = macroF1;
end