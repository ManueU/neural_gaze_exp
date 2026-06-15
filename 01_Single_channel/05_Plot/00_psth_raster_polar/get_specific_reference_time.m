function t_ref = get_specific_reference_time(data, bin_size, state_name)
%GET_SPECIFIC_REFERENCE_TIME Return onset time for one task state.
% If the state is missing, returns 0 and warns.

    [state_names, state_onset_s] = get_task_state_onsets(data, bin_size);
    idx = find(strcmpi(state_names, state_name), 1, 'first');

    if isempty(idx)
        t_ref = 0;
        warning('State "%s" not found. Using t_ref = 0.', state_name);
    else
        t_ref = state_onset_s(idx);
    end
end
