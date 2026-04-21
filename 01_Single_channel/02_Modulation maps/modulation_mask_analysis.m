% =========================================================
% DESCRIZIONE:
% Questo script:
% 1) decide se il canale è modulante (almeno un target)
% 2) costruisce modulation_mask {1 x n_arrays} (1 = modulante, NaN = non)
% 3) plotta gli array in bianco/grigio con numero di canale
%
% Richiede nello workspace:
%   data, mean_baseline, std_baseline,
%   PRE, POST, bin_size, period_pre, period_post,
%   n_sets, n_arrays, n_channels, n_targets
% =========================================================


% clearvars -except mean_baseline_common std_baseline_common
% close all
clc

% filename = '../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat'; 
load(filename);


%% Parametri
sets_plot = [1,2,3,4,5,6];
n_sets = numel(sets_plot);

n_arrays = 2;
n_channels = 96; 
n_targets = 8;
n_trials = 16;
bin_size = 0.02;
period_pre = 0.1;
period_post = 0.5;

% PRE = "Pres12";
% POST = "Reach";


%% Stabilisco se l'elettrodo è modulante o meno sfruttando il segnale zscored
% filename = '../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat'; 
% load(filename);
% PRE = "Pres12";
% POST = "Reach";
% 
% pre_bins  = max(1, round(period_pre/bin_size));
% post_bins = max(1, round(period_post/bin_size));
% len_bins  = pre_bins + post_bins;
% 
% idx_pres  = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == PRE);
% idx_reach = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == POST);
% 
% if isempty(idx_pres) || isempty(idx_reach)
%     error('Stati PRE/POST non trovati: controlla PRE/POST');
% end
% 
% modulation = false(n_arrays*n_channels, n_targets);
% z_thresh = 1; % soglia zscore
% run_len  = 3; % numero di bin consecutivi
% 
% for array = 1:n_arrays
%     for channel = 1:n_channels
%         ch_global = (array-1)*n_channels + channel;
% 
%         for target = 1:n_targets
%             sum_fr   = zeros(len_bins,1);
%             n_trials = 0;
% 
%             for set_ = 1:n_sets
%                 set = sets_plot(set_);
%                 idx = find([data(set).Data(array).Interp.Target_ID] == target & [data(set).Data(1).Interp.Excluded] == 0); 
%                 for k = 1:length(idx)
%                     tmp_pre = data(set).Data(array).Interp(idx(k)).Task_states{idx_pres, 2}(end-pre_bins+1:end, channel); 
%                     tmp_post = data(set).Data(array).Interp(idx(k)).Task_states{idx_reach,2}(1:post_bins, channel); 
%                     trial = [tmp_pre; tmp_post];  
%                     firing_rate = trial ./ bin_size;
% 
%                     sum_fr   = sum_fr + firing_rate;
%                     n_trials = n_trials + 1;
%                 end
%             end
% 
%             mean_across_trials = sum_fr / n_trials;   
% 
%             % z-score rispetto a baseline 
%             data_z = (mean_across_trials - mean_baseline_common(ch_global)) ./ std_baseline_common(ch_global);
% 
%             % maschera bin "significativi" 
%             mask = abs(data_z) > z_thresh;
% 
%             % almeno 3 bin consecutivi significativi?
%             modulation(ch_global, target) = any(strfind(mask', ones(1,run_len)));
% 
%         end
%     end
% end
% 
% % canali modulanti su almeno un target
% has_modulation = any(modulation, 2);  


%% Stabilisco se l'elettrodo è modulante o meno sfruttando PSTH su free-gaze
% load('responsive_channels_free_gaze.mat');
load(file)
responsive_channels = unique(results.channel_global); 
has_modulation = zeros(n_arrays*n_channels);
has_modulation(responsive_channels) = 1; 


%% Carico la mappa canali e costruisco modulation_mask
load('ChannelMap_BCI02.mat'); 

motor_medial   = ChannelMap.ChannelNumbers{1,1};
motor_lateral  = ChannelMap.ChannelNumbers{1,3};
motor_electrodes = {motor_medial, motor_lateral};  

modulation_mask = cell(1, n_arrays);

for array = 1:n_arrays
    M = motor_electrodes{array};   

    [n_rows, n_cols] = size(M);
    mask_array = nan(n_rows, n_cols);

    for i = 1:n_rows
        for j = 1:n_cols
            channel = M(i,j);   

            if isnan(channel)
                continue;
            end

            if has_modulation(channel)
                mask_array(i,j) = 1;     % modulante
            else
                mask_array(i,j) = NaN;   % non modulante
            end
        end
    end

    modulation_mask{array} = mask_array;
end


% %% Figure
% for array = 1:n_arrays
% 
%     M_elec = motor_electrodes{array};
%     M_mask = modulation_mask{array};
% 
%     [n_rows, n_cols] = size(M_elec);
% 
%     figure('Color','w'); hold on;
%     axis equal;
%     axis ij;  
% 
%     xlim([0.5, n_cols + 0.5]);
%     ylim([0.5, n_rows + 0.5]);
% 
%     clear set
%     set(gca, 'XTick', 1:n_cols, 'YTick', 1:n_rows);
%     box on;
% 
%     title(sprintf('Array %d', array));
% 
%     for i = 1:n_rows
%         for j = 1:n_cols
% 
%             ch = M_elec(i,j);
% 
%             if isnan(ch)
%                 continue;
%             end
% 
%             if ~isnan(M_mask(i,j)) && M_mask(i,j) == 1
%                 faceColor = [1 1 1];       % bianco se modulante
%             else
%                 faceColor = [0.7 0.7 0.7]; % grigio se non modulante
%             end
% 
%             rectangle('Position',[j-0.5, i-0.5, 1, 1], ...
%                       'FaceColor', faceColor, ...
%                       'EdgeColor', 'k');
% 
%             text(j, i, num2str(ch), ...
%                  'HorizontalAlignment','center', ...
%                  'VerticalAlignment','middle', ...
%                  'FontSize',8);
%         end
%     end
% end
% 
% disp('build_modulation_mask: calcolo e plot completati.');

