function response = find_response_peaks(t_win, y_win, baseline_val, thr, min_duration_bins)

    % inizializzazione output
    response = struct();

    response.baseline = baseline_val;
    response.threshold = thr;

    response.show_exc = false;
    response.show_inh = false;

    response.peak_exc_abs = NaN;
    response.peak_exc_rel = NaN;
    response.peak_exc_time = NaN;

    response.peak_inh_abs = NaN;
    response.peak_inh_rel = NaN;
    response.peak_inh_time = NaN;

    response.valid_exc = [];
    response.valid_inh = [];

    response.response_type = "none";

    % --------------------------------------------------------
    % Maschere sopra/sotto soglia
    % --------------------------------------------------------
    exc_mask = y_win >= (baseline_val + thr);
    inh_mask = y_win <= (baseline_val - thr);

    % segmenti consecutivi
    exc_runs = local_find_runs_internal(exc_mask);
    inh_runs = local_find_runs_internal(inh_mask);

    % --------------------------------------------------------
    % Tieni solo segmenti abbastanza lunghi
    % --------------------------------------------------------
    valid_exc = [];
    for rr = 1:size(exc_runs,1)
        if (exc_runs(rr,2) - exc_runs(rr,1) + 1) >= min_duration_bins
            valid_exc = [valid_exc; exc_runs(rr,:)]; 
        end
    end

    valid_inh = [];
    for rr = 1:size(inh_runs,1)
        if (inh_runs(rr,2) - inh_runs(rr,1) + 1) >= min_duration_bins
            valid_inh = [valid_inh; inh_runs(rr,:)]; 
        end
    end

    response.valid_exc = valid_exc;
    response.valid_inh = valid_inh;

    response.show_exc = ~isempty(valid_exc);
    response.show_inh = ~isempty(valid_inh);

    % --------------------------------------------------------
    % Picco eccitatorio
    % --------------------------------------------------------
    if response.show_exc
        best_exc_abs = -Inf;
        best_exc_time = NaN;
        best_exc_rel = NaN;

        for rr = 1:size(valid_exc,1)
            idx1 = valid_exc(rr,1);
            idx2 = valid_exc(rr,2);

            [tmp_abs, tmp_idx_local] = max(y_win(idx1:idx2));
            tmp_idx = idx1 + tmp_idx_local - 1;
            tmp_time = t_win(tmp_idx);
            tmp_rel = tmp_abs - baseline_val;

            if tmp_abs > best_exc_abs
                best_exc_abs = tmp_abs;
                best_exc_time = tmp_time;
                best_exc_rel = tmp_rel;
            end
        end

        response.peak_exc_abs = best_exc_abs;
        response.peak_exc_rel = best_exc_rel;
        response.peak_exc_time = best_exc_time;
    end

    % --------------------------------------------------------
    % Picco inibitorio
    % --------------------------------------------------------
    if response.show_inh
        best_inh_abs = Inf;
        best_inh_time = NaN;
        best_inh_rel = NaN;

        for rr = 1:size(valid_inh,1)
            idx1 = valid_inh(rr,1);
            idx2 = valid_inh(rr,2);

            [tmp_abs, tmp_idx_local] = min(y_win(idx1:idx2));
            tmp_idx = idx1 + tmp_idx_local - 1;
            tmp_time = t_win(tmp_idx);
            tmp_rel = tmp_abs - baseline_val;

            if tmp_abs < best_inh_abs
                best_inh_abs = tmp_abs;
                best_inh_time = tmp_time;
                best_inh_rel = tmp_rel;
            end
        end

        response.peak_inh_abs = best_inh_abs;
        response.peak_inh_rel = best_inh_rel;
        response.peak_inh_time = best_inh_time;
    end

    % --------------------------------------------------------
    % Tipo risposta
    % --------------------------------------------------------
    if response.show_exc && response.show_inh
        response.response_type = "both";
    elseif response.show_exc
        response.response_type = "exc";
    elseif response.show_inh
        response.response_type = "inh";
    else
        response.response_type = "none";
    end

end

function runs = local_find_runs_internal(mask)
    mask = mask(:)';
    dmask = diff([false, mask, false]);
    starts = find(dmask == 1);
    ends   = find(dmask == -1) - 1;
    runs = [starts(:), ends(:)];
end