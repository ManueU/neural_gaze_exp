function activity = compute_window_activity(trial_mat, bin_size, analysis_win_s, baseline_hz)
    if isempty(trial_mat) || any(isnan(analysis_win_s))
        activity = NaN;
        return
    end

    t = (0:size(trial_mat,1)-1) * bin_size;
    idx_win = t >= analysis_win_s(1) & t < analysis_win_s(2);

    if any(idx_win)
        activity = mean(trial_mat(idx_win,:), 'all') / bin_size - baseline_hz;
    else
        activity = NaN;
    end
end