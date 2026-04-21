function [pref_angle_map, tuning_strength_map] = gradientTargetPreferenceCircular(smoothed_modulation_matrix, data, labels_cond, d, valid_outline)

    if data.array == 1
        array_name = "Medial";
    else
        array_name = "Lateral";
    end

    n_targets = data.n_targets;
    angles = linspace(0, 2*pi, n_targets+1);
    angles(end) = [];

    map_size = size(smoothed_modulation_matrix{1,1});
    pref_angle_map = nan(map_size);
    tuning_strength_map = nan(map_size);

    for i = 1:map_size(1)
        for j = 1:map_size(2)

            tmp_mod_values = nan(1, n_targets);
            for target = 1:n_targets
                tmp_mod_values(target) = smoothed_modulation_matrix{1,target}(i, j);
            end

            valid = ~isnan(tmp_mod_values);
            if any(valid)
                original_vals = tmp_mod_values(valid);
                vals = original_vals - min(original_vals);
                angs = angles(valid);

                vec = nansum(vals .* exp(1i * angs));
                pref_angle_map(i,j) = angle(vec);
                tuning_strength_map(i,j) = abs(vec);
            end
        end
    end

    
    
    % porta gli angoli in [0, 2pi]
    pref_angle_plot = pref_angle_map;
    pref_angle_plot(pref_angle_plot < 0) = pref_angle_plot(pref_angle_plot < 0) + 2*pi;

    figure('Color','w')
    imagesc(pref_angle_plot, 'AlphaData', ~isnan(pref_angle_plot));
    hold on
    draw_array_outline(valid_outline);
    hold off
    axis image off
    colormap(slanCM(144))
    cb = colorbar;
    clim([0 2*pi])

    cb.Ticks = linspace(0, 2*pi, n_targets+1);
    cb.TickLabels = {'1','2','3','4','5','6','7','8','1'};

    set(gca, 'Color', [0.7 0.7 0.7])
    title({['\bf' char(labels_cond{d})], ['\rmArray ' char(array_name)]}, 'Interpreter','tex');
end