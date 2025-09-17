close all
clc

%% Targets
targ_coordinates = [Data.TaskStateMasks.uni_targets(1:9,2), Data.TaskStateMasks.uni_targets(1:9,3), Data.TaskStateMasks.uni_targets(1:9,1)];

x = targ_coordinates(:,1); 
y = targ_coordinates(:,2); 
figure('Color','w')
plot(x,y,'.', "MarkerSize", 10, "Color", "r")
axis equal
grid on
xlabel("x")
ylabel("y")
title('Targets position','Units','normalized','Position',[0.5 1.05 0]);

for i = 1:numel(x)
    text(x(i), y(i), sprintf('%d', i), ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','right', ...
        'FontSize',8, 'Color','k');
end


%% Single channel - timeXtrials (target)
ch = 40; 
n_targets = 8; 
bin_size = 0.02; 

for i = 1:n_targets
    M_spikes_tmp = [];
    for set = 1:4
        idx = find([dataset(set).Data(2).Resampled.Target_ID] == i);
        for j = 1:length(idx)
            M_spikes = [M_spikes_tmp, [dataset(set).Data(2).Resampled(idx(j)).Trial(:,ch)]]; 
            M_spikes_tmp = M_spikes;
        end 
    end 
    M_spikes_mean = mean(M_spikes, 2); 
    M_spikes_std  = std(M_spikes, 0, 2); % provare sem: mean standard error 

    firing_rate = M_spikes_mean ./ bin_size;
    firing_std  = M_spikes_std  ./ bin_size;  

    win_sec = 0.2; 
    w = max(3, round(win_sec / bin_size));   % numero di punti nella finestra
    fr_s   = smoothdata(firing_rate, 'gaussian', w);
    std_s  = smoothdata(firing_std,  'gaussian', w);

    t = (1:numel(firing_rate)) * bin_size;

    figure('Color','w'); hold on
    upper = fr_s + std_s;
    lower = fr_s - std_s;

    fill([t fliplr(t)], [upper' fliplr(lower')], ...
        [0.8 0.8 1], 'EdgeColor','none', 'FaceAlpha',0.3);

    plot(t, fr_s, 'b', 'LineWidth', 1.5)
    xline((29+55)*bin_size); 
    xline((29+55+54)*bin_size); 
    xline((29+55+54+156)*bin_size); 
    xline((29+55+54+156+4)*bin_size); 
    xline((29+55+54+156+4+32)*bin_size); 

    xlabel("Time (s)")
    ylabel("Firing rate (Hz)")
    title(['Channel = ' num2str(ch) '; Target = ' num2str(i)])
    hold off
end 
