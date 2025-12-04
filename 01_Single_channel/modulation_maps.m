clearvars -except mean_baseline std_baseline
% close all
clc

%% Modulation maps

filename = '../00_Data_extraction/free-gaze_BCI02.mat'; 
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
POST = "Reach";

%% Capire se l'elettrodo è modulante 
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
        mean_ch_tgt_mask = [];  
        trial_count = 0; 
        for target = 1:n_targets
            matrix_ch_tgt = [];
            for set = 1:n_sets
                idx = find([data(set).Data(array).Interp.Target_ID] == target); 
                for j = 1:length(idx)
                    trial_count = trial_count + 1; 
                    tmp_pre = data(set).Data(array).Interp(idx(j)).Task_states{idx_pres, 2}(end-pre_bins+1:end, channel); 
                    tmp_post = data(set).Data(array).Interp(idx(j)).Task_states{idx_reach,2}(1:post_bins, channel); 
                    trial = [tmp_pre; tmp_post];  

                    firing_rate = trial ./ bin_size;
                    
                    % firing_rate = trial ./ bin_size;
                    % trial_zscored = (firing_rate-mean_baseline(channel, array))./std_baseline(channel, array); 
                    % [~,id] = max(abs(trial_zscored));
                    % trial_max(trial_count, 1) = trial_zscored(id); 
                    % trial_label(trial_count, 1) = target;  

                    matrix_ch_tgt = [matrix_ch_tgt, firing_rate]; 
                    matrix_ch_all = [matrix_ch_all, firing_rate]; 
                end
            end
            mean_ch_tgt_mask(:,target, channel) = mean(matrix_ch_tgt, 2); 
            mean_ch_tgt(target, (array-1)*n_channels + channel) = mean(mean(matrix_ch_tgt));

        end 
        mean_all((array-1)*n_channels + channel) = mean(mean(matrix_ch_all, 2));
        % p_value((array-1)*n_channels + channel, 1) = anova1(mean_ch_tgt_mask_bs, label_mask_bs, 'off'); 
    end 
end 



%% Creazione di modulation matrix e modulation mask
load('ChannelMap_BCI02.mat'); 
motor_medial = ChannelMap.ChannelNumbers{1,1}; 
motor_lateral = ChannelMap.ChannelNumbers{1,3};
motor_electrodes = {motor_medial, motor_lateral};  

modulation_matrix = {nan(10,10), nan(10,10)}; 
modulation_mask = {nan(10,10), nan(10,10)}; 

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

            % Modulation mask
            if p_value(motor_electrodes{1,array}(i,j)) < 0.05
                modulation_mask{1,array}(i,j) = true;  
            else 
                modulation_mask{1,array}(i,j) = nan;  
            end 
            
        end 
    end 
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


