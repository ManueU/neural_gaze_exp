%% mappa somatotopica inspired by Natalya's paper
function [spearman_r] = gradientTargetPreference(smoothed_modulation_matrix, data, labels_cond, d)

    % --- colorbar ---
    start_color = [19, 31, 88]/255;    
    mid_color   = [87, 132, 165]/255;  
    end_color   = [148, 201, 187]/255;
    
    n_colors = 256;
    half = floor(n_colors/2);
    
    cmap1 = zeros(half,3);
    for i = 1:3
        cmap1(:,i) = linspace(start_color(i), mid_color(i), half);
    end
    
    cmap2 = zeros(n_colors - half, 3);
    for i = 1:3
        cmap2(:,i) = linspace(mid_color(i), end_color(i), n_colors - half);
    end
    
    cmap = [cmap1; cmap2];
    
    % --- parameters --- 
    if data.array == 1
        array = "Medial"; 
    else
        array = "Lateral";
    end 
    n_targets = data.n_targets;    
    target_identity = 1:n_targets;  
    
    map_size = size(smoothed_modulation_matrix{1,1});
    spearman_r = nan(map_size); 
    
    % --- compute Spearman correlation --- 
    for i = 1:map_size(1)
        for j = 1:map_size(2)
            tmp_mod_values = nan(1, n_targets);
            for target = 1:n_targets
                tmp_mod_values(target) = smoothed_modulation_matrix{1,target}(i, j);
            end
            % calcolo Spearman solo se non ci sono NaN
            if ~all(isnan(tmp_mod_values))
                spearman_r(i,j) = corr(tmp_mod_values(:), target_identity(:), 'Type', 'Spearman');
            end
        end
    end
    
    % ---- FIGURE POSITION ----
    figure('Color','w')
    
    h = imagesc(spearman_r, 'AlphaData', ~isnan(spearman_r)); 
    % colormap(cmap);
    colormap(slanCM(144))
    colorbar;
    axis equal tight;
    axis image 
    set(gca, 'XTick', [], 'YTick', [])
    nanColor = [0.7 0.7 0.7];   % grigio
    set(gca, 'Color', nanColor);
    clim([-1 1])
    
    title({['\bf' char(labels_cond{d})], ['\rmArray ' char(array)]}, 'Interpreter','tex');
end 


