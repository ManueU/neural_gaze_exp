function trial_cache = build_trial_cache(cond, target_ids, sets_plot, array_sel, channel_sel, bin_size)
    empty_cache = struct('trial_mat', [], 'raster_x', [], 'raster_y', []);
    trial_cache = repmat({empty_cache}, numel(cond), numel(target_ids));

    for c = 1:numel(cond)
        for tt = 1:numel(target_ids)
            [trial_mat, raster_x, raster_y] = collect_trials_for_target( ...
                cond(c).data, sets_plot, array_sel, channel_sel, target_ids(tt), bin_size);

            trial_cache{c,tt} = struct( ...
                'trial_mat', trial_mat, ...
                'raster_x', raster_x, ...
                'raster_y', raster_y);
        end
    end
end