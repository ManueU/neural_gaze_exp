close all
clc

filename = "../00_Data_extraction/controlled_BCI02.mat";
load(filename);

%% PARAMETRI
bin_size = 0.02;
n_sets   = 6;
channel  = 70;
array    = 2;      

target = 3; 
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

%% CALCOLO increment_times e labels 
% Se li hai già in workspace, puoi commentare questo blocco.
events_time_tmp = []; 
for i = 1:length(data(1).Data(array).Interp(1).Task_states)
    events_time = [events_time_tmp; ...
        size(data(1).Data(array).Interp(1).Task_states{i,2},1) * bin_size];
    events_time_tmp = events_time; 
end 
increment_times = cumsum(events_time); 
labels = string(data(1).Data(array).Interp(1).Task_states(:,1));

%% RACCOLTA TRIAL DEL SOLO TARGET
M_spikes   = [];
Trials_list = {};

for set = 1:n_sets
    idx = find([data(set).Data(array).Interp.Target_ID] == target);
    for j = 1:length(idx)
        trial_spikes = data(set).Data(array).Interp(idx(j)).Trial(:, channel);
        M_spikes = [M_spikes, trial_spikes];
        Trials_list{end+1} = trial_spikes;   %#ok<AGROW>
    end
end

%% CALCOLO PSTH (firing rate medio + deviazione standard)
M_mean   = mean(M_spikes, 2);
M_std    = std(M_spikes, 0, 2);
n_trials = size(M_spikes, 2);

firing_rate = M_mean ./ bin_size;   % spk/s
firing_std  = M_std  ./ bin_size;   % std del firing rate

% "Normalized" firing rate (sottraggo la baseline, in spk/s)
data_z = firing_rate - mean_baseline(channel, array);

% smoothing
w = 25;
data_z_s   = smoothdata(data_z,     'gaussian', w);  % media smussata
data_std_s = smoothdata(firing_std, 'gaussian', w);  % std smussata

% bande superiore/inferiore (± 1 std)
upper = data_z_s + data_std_s;
lower = data_z_s - data_std_s;

% asse dei tempi
n_bins = length(firing_rate);
t      = (1:n_bins) * bin_size;

%% FIGURA: 1) PSTH sopra (con banda std), 2) Raster sotto
figure('Color','w');
tiledlayout(2,1,"TileSpacing","compact","Padding","compact");

% ====================================================
% 1) PSTH NORMALIZZATA + DEVIAZIONE STANDARD
% ====================================================
ax1 = nexttile; hold on;

col = colors_target(target,:);

% banda ± std (semistrasparente)
tt      = t(:)';          % sicuro che sia riga
upperR  = upper(:)';
lowerR  = lower(:)';

fill([tt fliplr(tt)], [upperR fliplr(lowerR)], col, ...
     'FaceAlpha', 0.2, ...
     'EdgeColor', 'none', ...
     'HandleVisibility', 'off');   % non in legenda

% curva della media
plot(t, data_z_s, 'LineWidth', 1.8, ...
     'Color', col, ...
     'DisplayName', sprintf('Target %d', target));

if exist('increment_times','var')
    xline(increment_times,'k','HandleVisibility','off');
end

ylabel('Normalized Firing Rate');
title(sprintf('PSTH (Array %s, Ch %d) — Target %d', ...
      array_names(array), channel, target));
box on;
legend('Location','best');

% etichette dei task states
if exist('increment_times','var') && exist('labels','var')
    ax = gca;
    x_times = [0; increment_times(:)];
    x_text  = x_times(1:end-1) + diff(x_times)/2;
    y_text  = (ax.YLim(2) - 0.2) * ones(size(x_text));
    text(x_text, y_text, labels, 'HorizontalAlignment','center');
end

% ====================================================
% 2) RASTER SOTTO
% ====================================================
ax2 = nexttile; hold on;

t_edges    = (0:n_bins) * bin_size;
trial_counter = 0;
this_color = col;

for tr = 1:n_trials
    trial_counter = trial_counter + 1;
    trial_spikes  = Trials_list{tr};

    for b = 1:n_bins
        spike_count = trial_spikes(b);
        if spike_count == 0
            continue;
        end

        % distribuisco gli spike nel bin
        spike_times = linspace(t_edges(b), t_edges(b+1), spike_count + 2);
        spike_times = spike_times(2:end-1);

        % lineette verticali
        for st = spike_times
            plot([st st], [trial_counter-0.4 trial_counter+0.4], ...
                 'Color', this_color);
        end
    end
end

xlabel('Time (s)');
ylabel('Trial');
title('Raster');
clear set
set(gca,'YDir','reverse');
box off;

% sincronizza gli assi del tempo
linkaxes([ax1 ax2],'x');
