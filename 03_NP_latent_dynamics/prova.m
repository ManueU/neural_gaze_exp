nPC = 5; % usa più di 3 PC!

alpha_all = zeros(nCond-1, n_targets);
R2_all = zeros(nCond-1, n_targets);

for c = 2:nCond   % confronto con la prima condizione (baseline)
    for t = 1:n_targets
        
        nbin = size(condition_sep{1},1) / n_targets;
        
        rows1 = offset(1) + (t-1)*nbin + (1:nbin);
        rows2 = offset(c) + (t-1)*nbin + (1:nbin);
        
        traj1 = score(rows1,1:nPC);
        traj2 = score(rows2,1:nPC);

        % flatten (tempo × PC → vettore)
        x = traj1(:);
        y = traj2(:);

        % regressione senza intercetta
        alpha = (x' * y) / (x' * x);

        y_pred = alpha * x;
        R2 = 1 - sum((y - y_pred).^2) / sum((y - mean(y)).^2);

        alpha_all(c-1,t) = alpha;
        R2_all(c-1,t) = R2;
    end
end

disp('Scaling factors (alpha):')
disp(alpha_all)

disp('R^2 (quanto è puro scaling):')
disp(R2_all)

%%
figure; hold on

for c = 1:nCond
    for t = 1:n_targets
        
        nbin = size(condition_sep{c},1) / n_targets;
        rows = offset(c) + (t-1)*nbin + (1:nbin);
        
        traj = score(rows,1:5); % prime 5 PC
        norm_traj = vecnorm(traj,2,2);
        
        plot(norm_traj, 'LineWidth',1.5)
    end
end

title('Norma delle traiettorie PCA nel tempo')
xlabel('Time bins')
ylabel('||score||')


%%
for c = 2:nCond
    for t = 1:n_targets
        
        nbin = size(condition_sep{1},1) / n_targets;
        
        rows1 = offset(1) + (t-1)*nbin + (1:nbin);
        rows2 = offset(c) + (t-1)*nbin + (1:nbin);
        
        traj1 = score(rows1,1:5);
        traj2 = score(rows2,1:5);

        corr_val = corr(traj1(:), traj2(:));
        
        fprintf('Cond %d Target %d: corr = %.2f\n', c, t, corr_val);
    end
end

%%
figure; hold on

for c = 1:nCond
    mean_norm = zeros(nbin,1);
    
    for t = 1:n_targets
        rows = offset(c) + (t-1)*nbin + (1:nbin);
        traj = score(rows,1:5);
        norm_traj = vecnorm(traj,2,2);
        
        mean_norm = mean_norm + norm_traj;
    end
    
    mean_norm = mean_norm / n_targets;
    
    plot(mean_norm, 'LineWidth',3)
end

legend(labels)
title('Norma media delle traiettorie per condizione')