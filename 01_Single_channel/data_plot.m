% close all
clear all
clc

load("free-gaze_BCI02.mat");

%% Single channel
n_sets = 6; 
target_des = [1]; 
ch_start = 72;
ch_end = 72;
bin_size = 0.02;

events_time_tmp = []; 
for i = 1:length(data(1).Data(2).Resampled(1).Task_states)
    events_time = [events_time_tmp; size(data(1).Data(2).Resampled(1).Task_states{i,2},1)*bin_size];
    events_time_tmp = events_time; 
end 
increment_times = cumsum(events_time); 

labels = string(data(1).Data(2).Resampled(1).Task_states(:,1));

array_names = ["medial", "lateral"]; 
colors_target = [
    0.839, 0.153, 0.157;  % rosso
    0.122, 0.467, 0.706;  % blu
    0.172, 0.627, 0.172;  % verde
    0.580, 0.404, 0.741;  % viola
    1.000, 0.498, 0.055;  % arancione
    0.737, 0.741, 0.133;  % giallo oliva
    0.549, 0.337, 0.294;  % marrone
    0.890, 0.466, 0.760;  % rosa
];
% for array = 1:2
%     for channel = ch_start:ch_end
%         flag = 0; 
% 
%         for target = 1:length(target_des)
%             M_spikes_tmp = [];
% 
%             % M_spikes contains the spike count matrix of each trial
%             % associated to the specific target id. We consider all the
%             % sets.
%             for set = 1:n_sets
%                 idx = find([data(set).Data(array).Resampled.Target_ID] == target_des(target));                      
%                 for j = 1:length(idx)
%                     M_spikes = [M_spikes_tmp, [data(set).Data(array).Resampled(idx(j)).Trial(:,channel)]];   
%                     M_spikes_tmp = M_spikes;
%                 end
%             end 
%             M_spikes_mean = mean(M_spikes, 2); 
%             M_spikes_std  = std(M_spikes_mean);
%             M_spikes_sem  = std(M_spikes_mean)/sqrt(length(M_spikes_mean));
% 
%             firing_rate = M_spikes_mean ./ bin_size;
%             firing_std  = M_spikes_std  ./ bin_size;  
%             firing_sem  = M_spikes_sem  ./ bin_size;  
% 
%             mean_firing_rate = mean(firing_rate);
% 
%             w = 25; 
%             fr_s   = smoothdata(firing_rate, 'gaussian', w);
%             if(max(fr_s) > 2)
%                 if flag == 0
%                     figure('Color','w'); hold on
%                     if exist('increment_times','var') && ~isempty(increment_times)
%                       xline(increment_times, 'k', 'HandleVisibility','off');
%                     end
%                     yline(mean_firing_rate, 'Linewidth', 2, 'Color', [0.7, 0.7, 0.7], 'HandleVisibility','off');
% 
%                     xlabel('Time (s)');
%                     ylabel('Firing rate (Hz)');
%                     title(sprintf('Array = %s; Channel = %d;', array_names(array), channel));
%                     legend('Location','best');
%                     flag = 1; 
%                 end 
%                 t = (1:numel(firing_rate)) * bin_size;
%                 upper = fr_s + firing_std;
%                 lower = fr_s - firing_std;
% 
%                 fill([t fliplr(t)], [upper' fliplr(lower')], ...
%                     colors_target(target,:), 'EdgeColor','none', 'FaceAlpha',0.3, 'HandleVisibility','off');
%                 plot(t, fr_s, 'b', 'LineWidth', 1.5, 'Color', colors_target(target,:), 'DisplayName', sprintf('Target %d', target_des(target))), hold on
% 
%             end 
%         end 
%         if flag == 1
%             ax = gca;
%             x_times = [0; increment_times];
%             x_text = x_times(1:end-1) + diff(x_times)/2;
% 
%             y_text =  (ax.YLim(2) - 2)*ones(1,length(x_text)); 
%             text(x_text, y_text, labels, 'HorizontalAlignment', 'center'); 
%         end
%     end 
% end 


for array = 2:2
    for channel = ch_start:ch_end
        flag = 0; 
        
        for target = 1:length(target_des)
            M_spikes_tmp = [];

            % M_spikes contains the spike count matrix of each trial
            % associated to the specific target id. We consider all the
            % sets.
            for set = 1:n_sets
                idx = find([data(set).Data(array).Resampled.Target_ID] == target_des(target));                      
                for j = 1:length(idx)
                    M_spikes = [M_spikes_tmp, [data(set).Data(array).Resampled(idx(j)).Trial(:,channel)]];   
                    M_spikes_tmp = M_spikes;
                end
            end 
            M_spikes_mean = mean(M_spikes, 2); 
            M_spikes_std  = std(M_spikes_mean);
            M_spikes_sem  = std(M_spikes_mean)/sqrt(length(M_spikes_mean));
        
            firing_rate = M_spikes_mean ./ bin_size;
            firing_std  = M_spikes_std  ./ bin_size;  
            firing_sem  = M_spikes_sem  ./ bin_size;  

            mean_firing_rate = mean(firing_rate);
            data_zscored = (firing_rate-mean_firing_rate)./firing_std; 
            
            w = 25; 
            data_zscored_s = smoothdata(data_zscored, 'gaussian', w); 
            fr_s   = smoothdata(firing_rate, 'gaussian', w);
            if(max(fr_s) > 2)
                if flag == 0
                    figure('Color','w'); hold on
                    if exist('increment_times','var') && ~isempty(increment_times)
                      xline(increment_times, 'k', 'HandleVisibility','off');
                    end

                    xlabel('Time (s)');
                    ylabel('z-score');
                    title(sprintf('Array = %s; Channel = %d;', array_names(array), channel));
                    legend('Location','best');
                    flag = 1; 
                end 
                t = (1:numel(firing_rate)) * bin_size;
                plot(t, data_zscored_s, 'b', 'LineWidth', 1.5, 'Color', colors_target(target,:), 'DisplayName', sprintf('Target %d', target_des(target))), hold on
            
            end 
        end 
        if flag == 1
            ax = gca;
            x_times = [0; increment_times];
            x_text = x_times(1:end-1) + diff(x_times)/2;
            
            y_text =  (ax.YLim(2) - 0.2)*ones(1,length(x_text)); 
            text(x_text, y_text, labels, 'HorizontalAlignment', 'center'); 
        end
    end 
end 
