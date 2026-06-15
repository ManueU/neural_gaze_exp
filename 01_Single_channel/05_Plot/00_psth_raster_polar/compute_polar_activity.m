function polar_activity = compute_polar_activity(trial_cache, cond, cond_idx, target_ids, ...
        bin_size, win_rel, baseline_hz, ref_state)

    polar_activity = nan(numel(cond_idx), numel(target_ids));

    for ii = 1:numel(cond_idx)
        c = cond_idx(ii);
        t_ref = get_specific_reference_time(cond(c).data, bin_size, ref_state);
        analysis_win_s = t_ref + win_rel;

        for tt = 1:numel(target_ids)
            polar_activity(ii,tt) = compute_window_activity(trial_cache{c,tt}.trial_mat, bin_size, analysis_win_s, baseline_hz);
        end
    end
end