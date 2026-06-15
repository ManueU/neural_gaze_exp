function [event_times, event_labels] = get_condition_events(data, bin_size, cond_code)
%GET_CONDITION_EVENTS Return event times and labels for one condition.

    [state_names, state_onset_s] = get_task_state_onsets(data, bin_size);

    switch lower(string(cond_code))
        case {"free-gaze", "motor"}
            event_names  = ["pres12", "reach"];
            event_labels = ["Target", "Go"];

        case "controlled"
            event_names  = ["pres12", "gaze", "reach"];
            event_labels = ["Target", "Saccade", "Go"];

        case "gaze-only"
            event_names  = ["pres12", "gaze"];
            event_labels = ["Target", "Saccade"];

        otherwise
            error('Unknown condition code: %s', cond_code);
    end

    event_times = nan(size(event_names));
    for k = 1:numel(event_names)
        idx = find(strcmpi(state_names, event_names(k)), 1, 'first');
        if ~isempty(idx)
            event_times(k) = state_onset_s(idx);
        end
    end

    valid = ~isnan(event_times);
    event_times = event_times(valid);
    event_labels = event_labels(valid);
end
