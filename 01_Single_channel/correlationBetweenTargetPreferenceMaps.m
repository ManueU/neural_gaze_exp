%% correlation between  two smoothed_modulation_matrix
function correlationBetweenTargetPreferenceMaps(correlation_maps, data, labels_cond)

    if data.array == 1
       array = "Medial"; 
    else
       array = "Lateral";
    end 
    nMaps = numel(correlation_maps);
    pairs = nchoosek(1:nMaps, 2);
    nPairs = size(pairs, 1);
    
    num_shuffle = 10000;
    
    real_corr = nan(nPairs, 1);
    p_value = nan(nPairs, 1);
    n_valid = nan(nPairs, 1);
    shuffle_corr = nan(num_shuffle, nPairs);
    
    % --- compute the correlation ---
    for p = 1:nPairs
        i = pairs(p, 1);
        j = pairs(p, 2);
    
        map1 = correlation_maps{i};
        map2 = correlation_maps{j};
    
        v1 = map1(:);
        v2 = map2(:);
    
        valid = ~isnan(v1) & ~isnan(v2);
        n_valid(p) = sum(valid);
    
        if n_valid(p) < 2
            real_corr(p) = NaN; 
            p_value(p) = NaN;
            continue
        end 
        
        mod1_vector = v1(valid);
        mod2_vector = v2(valid);
        real_corr(p) = corr(mod1_vector, mod2_vector); 
        
        % Shuffling
        shuffle_corr = zeros(num_shuffle, 1);
        L = numel(mod2_vector);
        for s = 1:num_shuffle
            idx = randperm(L);
            shuffled_mod2 = mod2_vector(idx);
            shuffle_corr(s) = corr(mod1_vector, shuffled_mod2);
        end
        p_value(p) = mean(shuffle_corr >= real_corr(p));
    
        % Plot e p-value
        figure('Color','white');
        histogram(shuffle_corr, 'DisplayName','shuffle');
        hold on
        xline(real_corr(p), 'r', 'LineWidth', 2, 'DisplayName','real');
        xlabel('Correlation');
        ylabel('Frequency');
        title({ ...
            ['\bf' char(labels_cond{i}) ' vs ' char(labels_cond{j})], ...
            ['\rmArray ' char(array)]}, ...
            'Interpreter','tex');        xlim([-1 1]);
        legend('show');
    
    end 
end 
