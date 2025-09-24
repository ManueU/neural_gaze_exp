clc 
close all

%% Tuning curve X target
target_des = [2, 4, 6]; 
ch_start = 1;
ch_end = 96;
n_sets = 4; 
bin_size = 0.02;

 
% 1: medial arm 
% 2: lateral hand 
phases_names = string(dataset(1).Data(1).Resampled(1).Task_states(:,1)); 
array_names = ["medial", "lateral"]; 

for array = 2:2
    for channel = ch_start:ch_end
        flag = 0; 
        for phase = 1:length(phases_names)
            firing_rate_mean = []; 
            firing_sem = []; 
            for target = 1:length(target_des)
                M_spikes_tmp = [];
                for set = 1:n_sets
                    idx = find([dataset(set).Data(array).Resampled.Target_ID] == target_des(target));                      
                    for j = 1:length(idx)
                        M_spikes = [M_spikes_tmp, [dataset(set).Data(array).Resampled(idx(j)).Task_states{phase,2}(:,channel)]];   
                        M_spikes_tmp = M_spikes;
                    end
                end 
                M_spikes_mean = mean(M_spikes, 2); 
                M_spikes_std  = std(M_spikes_mean);
                M_spikes_sem  = std(M_spikes_mean)/sqrt(length(M_spikes_mean));
            
                firing_rate = M_spikes_mean ./ bin_size;
                firing_std  = M_spikes_std  ./ bin_size;  
                
                tmp_sem  = M_spikes_sem  ./ bin_size; 
                firing_sem = [firing_sem; tmp_sem]; 
                
                tmp_mean = mean(firing_rate);
                firing_rate_mean = [firing_rate_mean; tmp_mean];
            end
            if max(firing_rate_mean > 0)
                if flag == 0
                    figure('Color','w'); hold on
                    xlim([0, 9])
                    flag = 1; 
                end 
                errorbar(target_des, firing_rate_mean, firing_sem, '-o','MarkerSize', 8, 'LineWidth', 1.2, 'DisplayName', phases_names{phase}), hold on
                legend('show','Location','best');
            end 
        end 
        yline(baseline_const(channel, array), '--', 'HandleVisibility','off');
        % text(8.5, baseline_const(channel, array)+0.5, 'Baseline')
        xlabel('Target');
        ylabel('Mean firing rate (Hz)');
        title(sprintf('Array = %s; Channel = %d;', array_names(array), channel));

    end 
end 


%% Tuning curve X phase
clearvars -except dataset baseline baseline_const

target_des = [2, 4, 6, 8]; 
ch_start = 29;
ch_end = 29;
n_sets = 4; 
bin_size = 0.02;
colors_target = [ ...
    0.8  0.4  0.4;  % rosso 
    0.6  0.8  0.6;  % verde 
    0.4  0.6  0.8;  % azzurro 
    0.8  0.8  0.5]; % giallo 
 
% 1: medial arm 
% 2: lateral hand 
phases_names = string(dataset(1).Data(1).Resampled(1).Task_states(:,1)); 
phase_array = 1:length(phases_names); 
array_names = ["medial", "lateral"]; 

for array = 2:2
    for channel = ch_start:ch_end
        figure('Color','w'); hold on
        for target = 1:length(target_des) 
            firing_rate_mean = []; 
            firing_sem = []; 
            for phase = 1:length(phases_names)
                M_spikes_tmp = [];
                for set = 1:n_sets
                    idx = find([dataset(set).Data(array).Resampled.Target_ID] == target_des(target)); 
                    for j = 1:length(idx)
                        M_spikes = [M_spikes_tmp, [dataset(set).Data(array).Resampled(idx(j)).Task_states{phase,2}(:,channel)]];   
                        M_spikes_tmp = M_spikes;
                    end
                end
                M_spikes_mean = mean(M_spikes, 2); 
                M_spikes_std  = std(M_spikes_mean);
                M_spikes_sem  = std(M_spikes_mean)/sqrt(length(M_spikes_mean));
            
                firing_rate = M_spikes_mean ./ bin_size;
                firing_std  = M_spikes_std  ./ bin_size;  
                
                tmp_sem  = M_spikes_sem  ./ bin_size; 
                firing_sem = [firing_sem; tmp_sem]; 
                
                tmp_mean = mean(firing_rate);
                firing_rate_mean = [firing_rate_mean; tmp_mean];
            end 
            errorbar(phase_array, firing_rate_mean, firing_sem, '-o','MarkerSize', 8, 'LineWidth', 1.2, 'Color', colors_target(target, :), 'DisplayName', string(target_des(target))), hold on

        end 
        lgd = legend('show','Location','best');
        title(lgd, 'Target IDs');

        yline(baseline_const(channel, array), '--', 'HandleVisibility','off');
        xticks(1:length(phases_names));
        xticklabels(phases_names);
        xlim([0.5 length(phases_names)+0.5]);     
        xtickangle(0);            
        xlabel('Phase');
        ylabel('Mean firing rate (Hz)');
        title(sprintf('Array = %s; Channel = %d;', array_names(array), channel));
    end 
end 
