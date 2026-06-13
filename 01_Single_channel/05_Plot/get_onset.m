function t = get_onset(data, bin_size, state_name)

    % cerca il primo trial valido
    found = false;
    for s = 1:numel(data)
        for a = 1:numel(data(s).Data)
            if isfield(data(s).Data(a),'Interp') && ~isempty(data(s).Data(a).Interp)
                TS = data(s).Data(a).Interp(1).Task_states;
                found = true;
                break
            end
        end
        if found, break; end
    end

    if ~found
        error('Task_states non trovati');
    end

    state_names = string(TS(:,1));
    state_dur_s = cellfun(@(x) size(x,1)*bin_size, TS(:,2));
    state_onset_s = [0; cumsum(state_dur_s(1:end-1))];

    idx = find(strcmpi(state_names, state_name), 1, 'first');

    if isempty(idx)
        t = NaN;
    else
        t = state_onset_s(idx);
    end
end