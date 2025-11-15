function data_resampled = interp_data(data, len_max)
        n_orig   = 1:size(data,1);                              
        n_target = linspace(1, size(data,1), len_max);          
        data_resampled = round(interp1(n_orig, data, n_target, 'linear'));
end 