function acc = lda_cv(X, Y, kFolds)
    % Una stima di accuratezza con 5-fold CV (senza ripetizioni)
    cvp = cvpartition(Y, 'KFold', kFolds);
    Mdl = fitcdiscr(X, Y, 'DiscrimType','pseudolinear', 'Prior','empirical');
    CVMdl = crossval(Mdl, 'CVPartition', cvp);
    loss = kfoldLoss(CVMdl, 'LossFun','ClassifError');
    acc = 1 - loss;
end 