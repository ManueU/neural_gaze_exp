function [state_names, state_onset_s] = get_task_state_onsets(data, bin_size)
%GET_TASK_STATE_ONSETS Return task-state names and onset times in seconds.

    TS = [];

    for s = 1:numel(data)
        if ~isfield(data(s), 'Data')
            continue
        end

        for a = 1:numel(data(s).Data)
            if isfield(data(s).Data(a), 'Interp') && ~isempty(data(s).Data(a).Interp) && ...
                    isfield(data(s).Data(a).Interp(1), 'Task_states')
                TS = data(s).Data(a).Interp(1).Task_states;
                break
            end
        end

        if ~isempty(TS)
            break
        end
    end

    if isempty(TS)
        error('Task_states not found.');
    end

    state_names = string(TS(:,1));
    state_dur_s = cellfun(@(x) size(x,1) * bin_size, TS(:,2));
    state_onset_s = [0; cumsum(state_dur_s(1:end-1))];
end
