function mask_move = convert_p_to_mask(anova_map)
    mask_move = anova_map;
    mask_move(mask_move>=0.05) = nan;
    mask_move(mask_move<0.05) = true;
end
