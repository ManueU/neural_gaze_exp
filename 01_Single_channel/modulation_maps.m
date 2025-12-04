clearvars -except mean_baseline std_baseline
% close all
clc

%% Modulation maps

filename = '../00_Data_extraction/gaze_BCI02.mat'; 
load(filename);

n_sets = 6;
n_arrays = 2;
n_channels = 96; 
n_targets = 8; 
n_trials = 32;
bin_size = 0.02;
period_pre = 0.1;
period_post = 0.5;
PRE = "Pres12";
POST = "Gaze";

%% Modulation matrix
idx_pres = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == PRE); 
idx_reach = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == POST); 
if isempty(idx_pres) || isempty(idx_reach)
    error('Stati PRE/POST non trovati: controlla PRE/POST');
end

pre_bins  = max(1, round(period_pre/bin_size));
post_bins = max(1, round(period_post/bin_size));
for array = 1:n_arrays
    for channel = 1:n_channels
        matrix_ch_all = []; 
        for target = 1:n_targets
            matrix_ch_tgt = [];
            for set = 1:n_sets
                idx = find([data(set).Data(array).Interp.Target_ID] == target); 
                for j = 1:length(idx)
                    tmp_pre = data(set).Data(array).Interp(idx(j)).Task_states{idx_pres, 2}(end-pre_bins+1:end, channel); 
                    tmp_post = data(set).Data(array).Interp(idx(j)).Task_states{idx_reach,2}(1:post_bins, channel); 
                    trial = [tmp_pre; tmp_post];  

                    firing_rate = trial ./ bin_size;

                    matrix_ch_tgt = [matrix_ch_tgt, firing_rate]; 
                    matrix_ch_all = [matrix_ch_all, firing_rate]; 
                end
            end
            mean_ch_tgt(target, (array-1)*n_channels + channel) = mean(mean(matrix_ch_tgt));

        end 
        mean_all((array-1)*n_channels + channel) = mean(mean(matrix_ch_all, 2));
    end 
end 

%% Creazione di modulation matrix
load('ChannelMap_BCI02.mat'); 
motor_medial = ChannelMap.ChannelNumbers{1,1}; 
motor_lateral = ChannelMap.ChannelNumbers{1,3};
motor_electrodes = {motor_medial, motor_lateral};  

modulation_matrix = {nan(10,10), nan(10,10)}; 

for array = 1:n_arrays
    for i = 1:10
        for j = 1:10
            if isnan(motor_electrodes{1,array}(i,j))
                continue;
            end 
            
            % Modulation matrix
            for target = 1:n_targets
                modulation_matrix{target,array}(i,j) = (mean_ch_tgt(target, motor_electrodes{1,array}(i,j)) - mean_all(1,motor_electrodes{1,array}(i,j)))./mean_all(1,motor_electrodes{1,array}(i,j));
            end 

        end 
    end 
end


%% Modulation mask
for_mask  = cell(n_arrays, n_channels, n_targets);

for array = 1:n_arrays
    for channel = 1:n_channels
        for target = 1:n_targets
            matrix_for_mask = [];
            for set = 1:n_sets
                idx = find([data(set).Data(array).Interp.Target_ID] == target); 
                for j = 1:length(idx)
                    tmp_pre  = data(set).Data(array).Interp(idx(j)).Task_states{idx_pres, 2}(end-pre_bins+1:end, channel); 
                    tmp_post = data(set).Data(array).Interp(idx(j)).Task_states{idx_reach,2}(1:post_bins, channel); 
                    trial    = [tmp_pre; tmp_post];  

                    firing_rate = trial ./ bin_size;
                    matrix_for_mask = [matrix_for_mask, firing_rate]; 
                end
            end
            for_mask{array, channel, target} = matrix_for_mask;
        end
    end
end



%% Permutation test
n_perm = 10000;             % numero di permutazioni
alpha  = 0.05;              % soglia di significatività

p_perm   = nan(n_arrays*n_channels, n_targets);
sig_mask = false(n_arrays*n_channels, n_targets);

for array = 1:n_arrays
    for channel = 1:n_channels
        baseline = mean_baseline(channel, array); 
      
        for target = 1:n_targets
            X = for_mask{array, channel, target}; 
            trial_mean = mean(X, 1);                
            d = trial_mean - baseline;            
            
            % permutation test one-sample
            p = permutation_test_onesample(d, n_perm);
            
            p_perm((array-1)*n_channels + channel, target)   = p;
            sig_mask((array-1)*n_channels + channel, target) = (p < alpha);
        end
        
    end
end

%% Modulation mask a livello di canale
modulation_mask = cell(1, n_arrays);

for array = 1:n_arrays
    mask_array = nan(10,10);  

    for i = 1:10
        for j = 1:10

            channel = motor_electrodes{array}(i,j);  % numero del canale

            if isnan(channel)
                continue;
            end

            % controllo se almeno 1 target è significativo
            if any(sig_mask(channel, :), 'all')
                mask_array(i,j) = 1;   % modulante
            else
                mask_array(i,j) = nan; % non modulante
            end

        end
    end

    modulation_mask{array} = mask_array;
end




%% Color map
blue = [106/255, 174/255, 233/255];
white = [1, 1, 1];
red = [209/255, 69/255, 128/255];
black = [0.7 0.7 0.7];
map = [linspace(blue(1), white(1), 100)', linspace(blue(2), white(2), 100)', linspace(blue(3), white(3), 100)';
        linspace(white(1), red(1), 100)', linspace(white(2), red(2), 100)', linspace(white(3), red(3), 100)'];
map = [black; map];



%%
% for mod =  [1 2 3 4]  % : length(modality_name)
    for array = 2%: length(arrays_order)
        figure('Color', 'w')
        for target = 1:n_targets
            masked_index = modulation_matrix{target, array} .* modulation_mask{1,array};
            % smoothed_masked = nangauss_smooth(masked_index, 4, 0.5);
            subplot(1, n_targets, target)
            h = imagesc(masked_index)
            colormap(map)
            clim([-1 1]);
        end
    end
% end


