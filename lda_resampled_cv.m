function acc_dist = lda_resampled_cv(X, Y, K, n_rep)
    acc_dist  = zeros(n_rep,1);
    for r = 1:n_rep
        cv = cvpartition(Y,'KFold',K,'Stratify',true);
        fold_acc = zeros(cv.NumTestSets,1);
        for k = 1:cv.NumTestSets
            idxTrain = training(cv, k);
            idxTest  = test(cv, k);            
            % mu = mean(X(idxTrain,:),1);  
            % sg = std(X(idxTrain,:),[],1);
            % sg(sg==0) = 1; % evita divisione per zero
            % Xtrain = (X(idxTrain,:) - mu) ./ sg;
            % Xtest = (X(idxTest,:) - mu) ./ sg;

            Mdl = fitcdiscr(Xtrain, Y(idxTrain), 'DiscrimType','pseudoLinear'); 
            Ypred = predict(Mdl, Xtest);
            fold_acc(k) = mean(Ypred == Y(idxTest));
        end
        acc_dist(r) = mean(fold_acc);
    end
end