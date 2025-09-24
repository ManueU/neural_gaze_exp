% close all
clc


%% Single channel
target_des = 6; 
ch_start = 29;
ch_end = 29;
n_sets = 4; 
bin_size = 0.02;

events_time_tmp = []; 
for i = 1:length(dataset(1).Data(2).Resampled(1).Task_states)
    events_time = [events_time_tmp; size(dataset(1).Data(2).Resampled(1).Task_states{i,2},1)*bin_size];
    events_time_tmp = events_time; 
end 
increment_times = cumsum(events_time); 

labels = string(dataset(1).Data(2).Resampled(1).Task_states(:,1));

% 1: medial arm 
% 2: lateral hand 
array_names = ["medial", "lateral"]; 
colors_target = [ ...   
    0.4  0.6  0.8;    
    0.8  0.4  0.4;  % rosso 
    0.6  0.8  0.6;  % verde 
    0.4  0.6  0.8;  % azzurro
    0.8  0.8  0.5]; % giallo 
    

for array = 1:2
    for channel = ch_start:ch_end
        flag = 0; 
    
        for target = 1:length(target_des)
            M_spikes_tmp = [];
            for set = 1:n_sets
                idx = find([dataset(set).Data(array).Resampled.Target_ID] == target_des(target));                      
                for j = 1:length(idx)
                    M_spikes = [M_spikes_tmp, [dataset(set).Data(array).Resampled(idx(j)).Trial(:,channel)]];   
                    M_spikes_tmp = M_spikes;
                end
            end 
            M_spikes_mean = mean(M_spikes, 2); 
            M_spikes_std  = std(M_spikes_mean);
            M_spikes_sem  = std(M_spikes_mean)/sqrt(length(M_spikes_mean));
        
            firing_rate = M_spikes_mean ./ bin_size;
            firing_std  = M_spikes_std  ./ bin_size;  
            firing_sem  = M_spikes_sem  ./ bin_size;  
            
            w = 25; 
            fr_s   = smoothdata(firing_rate, 'gaussian', w);
            
            if(max(fr_s) > 2)
                if flag == 0
                    figure('Color','w'); hold on
                    if exist('increment_times','var') && ~isempty(increment_times)
                      xline(increment_times, 'k', 'HandleVisibility','off');
                    end
        
                    if array == 2
                        ch = channel + 96; 
                    else 
                        ch = channel; 
                    end 
                    bas_smooth = smoothdata(baseline(ch).mean, 'gaussian', w);
            
                    t_bas = 0:bin_size:increment_times(end)-bin_size; 
                    upper = bas_smooth + baseline(ch).std;
                    lower = bas_smooth - baseline(ch).std;
                    fill([t_bas fliplr(t_bas)], [upper' fliplr(lower')], 'k', ...
                    'EdgeColor','none', 'FaceAlpha',0.1, 'HandleVisibility','off');
                    plot(t_bas, bas_smooth, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'HandleVisibility','off');
                    % yline(baseline_const(channel,array), 'HandleVisibility','off');

                    xlabel('Time (s)');
                    ylabel('Firing rate (Hz)');
                    title(sprintf('Array = %s; Channel = %d;', array_names(array), channel));
                    legend('Location','best');
                    flag = 1; 
                end 
                t = (1:numel(firing_rate)) * bin_size;
                upper = fr_s + firing_std;
                lower = fr_s - firing_std;
            
                fill([t fliplr(t)], [upper' fliplr(lower')], ...
                    colors_target(target,:), 'EdgeColor','none', 'FaceAlpha',0.3, 'HandleVisibility','off');
                plot(t, fr_s, 'b', 'LineWidth', 1.5, 'Color', colors_target(target,:), 'DisplayName', sprintf('Target %d', target_des(target))), hold on
            
            end 
        end 
        if flag == 1
            ax = gca;
            x_times = [0; increment_times];
            x_text = x_times(1:end-1) + diff(x_times)/2;
            
            y_text =  (ax.YLim(2) - 2)*ones(1,length(x_text)); 
            text(x_text, y_text, labels, 'HorizontalAlignment', 'center'); 
        end
    end 
end 
