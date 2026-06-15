function t = get_state_onset(state_names, state_onset_s, name)
    idx = find(strcmpi(state_names, name), 1, 'first');
    if isempty(idx)
        t = nan;
    else
        t = state_onset_s(idx);
    end
end