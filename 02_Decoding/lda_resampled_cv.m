function acc_dist = lda_resampled_cv(X, Y, K, n_rep)
    acc_dist  = zeros(n_rep,1);
    for r = 1:n_rep
        fold_acc = lda_cv(X, Y, K); 
        acc_dist(r) = mean(fold_acc);
    end
end