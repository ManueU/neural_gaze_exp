function [trial_mat, raster_x, raster_y] = collect_trials_for_target(data, sets_plot, array_sel, channel_sel, target_id, bin_size)
%COLLECT_TRIALS_FOR_TARGET Collect spike trains and raster points for one target.

    trial_cols = {};
    raster_x = [];
    raster_y = [];
    trial_id = 0;

    for s = sets_plot(:)'
        if s > numel(data) || ~isfield(data(s), 'Data') || array_sel > numel(data(s).Data)
            continue
        end

        if ~isfield(data(s).Data(array_sel), 'Interp') || isempty(data(s).Data(array_sel).Interp)
            continue
        end

        trials = data(s).Data(array_sel).Interp;
        keep = find([trials.Target_ID] == target_id & [trials.Excluded] == 0);

        for k = keep(:)'
            if channel_sel > size(trials(k).Trial, 2)
                continue
            end

            spk = trials(k).Trial(:, channel_sel);
            if isempty(spk)
                continue
            end

            trial_id = trial_id + 1;
            trial_cols{end+1} = spk(:); 

            spike_bins = find(spk > 0);
            raster_x = [raster_x; (spike_bins - 1) * bin_size]; 
            raster_y = [raster_y; trial_id * ones(numel(spike_bins), 1)]; 
        end
    end

    if isempty(trial_cols)
        trial_mat = [];
        return
    end

    min_len = min(cellfun(@numel, trial_cols));
    trial_mat = cell2mat(cellfun(@(x) x(1:min_len), trial_cols, 'UniformOutput', false));
end
