close all
clc

filename = "../00_Data_extraction/free-gaze_BCI02.mat";
load(filename);

%% Parametri
bin_size = 0.02;
n_sets = 6;
channel = 64;
array = 2;      
targets = [2 3 4];

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

%% Calcolo increment_times (eventi di stato lungo il trial)
events_time_tmp = []; 
for i = 1:length(data(1).Data(array).Interp(1).Task_states)
    events_time = [events_time_tmp; ...
        size(data(1).Data(array).Interp(1).Task_states{i,2},1) * bin_size];
    events_time_tmp = events_time; 
end 
increment_times = cumsum(events_time);   % tempi assoluti (da inizio trial)

%% Raccolta trial per ogni target
M_spikes_all = cell(numel(targets),1);
Trials_list_all = cell(numel(targets),1);

for tIdx = 1:numel(targets)
    tgt = targets(tIdx);

    M_spikes    = [];
    Trials_list = {};

    for set = 1:n_sets
        idx = find([data(set).Data(array).Interp.Target_ID] == tgt);
        for j = 1:length(idx)
            trial_spikes = data(set).Data(array).Interp(idx(j)).Trial(:, channel);
            M_spikes = [M_spikes, trial_spikes];
            Trials_list{end+1} = trial_spikes;   
        end
    end

    M_spikes_all{tIdx}    = M_spikes;
    Trials_list_all{tIdx} = Trials_list;
end

%% Calcolo PSTH (firing rate medio, std, sem) + asse tempi allineato al Reach

% riferimento per numero di bin
M_spikes_ref = M_spikes_all{1};
n_bins = size(M_spikes_ref, 1);

% indice stato Reach
idx_reach = find(string(data(1).Data(array).Interp(1).Task_states(:,1)) == "Reach");

% tempo di onset del Reach (inizio dello stato Reach rispetto a inizio trial)
if idx_reach == 1
    reach_onset = 0;
else
    reach_onset = increment_times(idx_reach-1);  % inizio dello stato Reach
end

% asse tempi allineato al Reach
t  = ((1:n_bins) * bin_size) - reach_onset;
t_edges = (0:n_bins) * bin_size   - reach_onset;

% allineo anche i marker degli stati allo stesso riferimento
increment_times_aligned = increment_times - reach_onset;

firing_rate_all = cell(numel(targets),1);
firing_sem_all  = cell(numel(targets),1);
data_z_all      = cell(numel(targets),1);
upper_all       = cell(numel(targets),1);
lower_all       = cell(numel(targets),1);

w = 25; % smoothing window

for tIdx = 1:numel(targets)
    M_spikes = M_spikes_all{tIdx};
    if isempty(M_spikes)
        warning('Nessun trial per target %d', targets(tIdx));
        continue;
    end

    M_mean  = mean(M_spikes, 2);
    M_std   = std(M_spikes, 0, 2);
    nTrials = size(M_spikes, 2);
    M_sem   = M_std ./ sqrt(nTrials);

    firing_rate = M_mean ./ bin_size;   % spk/s
    firing_std  = M_std  ./ bin_size;
    firing_sem  = M_sem  ./ bin_size;

    % "Normalized" firing rate (sottraggo la baseline, in spk/s)
    data_z = firing_rate - mean_baseline(channel, array);

    % smoothing
    data_z_s   = smoothdata(data_z, 'gaussian', w);
    data_std_s = smoothdata(firing_std, 'gaussian', w);
    data_sem_s = smoothdata(firing_sem, 'gaussian', w);

    upper = data_z_s + data_std_s;
    lower = data_z_s - data_std_s;

    firing_rate_all{tIdx} = data_z_s;
    firing_sem_all{tIdx} = data_std_s;
    data_z_all{tIdx} = data_z_s;
    upper_all{tIdx} = upper;
    lower_all{tIdx} = lower;
end

%% Figure (1): 1) PSTH (multi-target), 2) Raster multi-target
figure('Color','w');
tiledlayout(2,1,"TileSpacing","compact","Padding","compact");

% ================== PSTH ==================
ax1 = nexttile; hold on;

for tIdx = 1:numel(targets)
    if isempty(data_z_all{tIdx}), continue; end

    tgt = targets(tIdx);
    col = colors_target(tgt,:);

    tt = t(:)';          
    upperR = upper_all{tIdx}(:)';
    lowerR = lower_all{tIdx}(:)';

    % banda ± sem
    fill([tt fliplr(tt)], [upperR fliplr(lowerR)], col, ...
         'FaceAlpha', 0.15, ...
         'EdgeColor', 'none', ...
         'HandleVisibility', 'off');  

    % curva media
    plot(t, data_z_all{tIdx}, 'LineWidth', 1.8, ...
         'Color', col, ...
         'DisplayName', sprintf('Target %d', tgt));
end

xline(0,'k','HandleVisibility','off');
clear set
set(gca, 'XTick', []);
ylabel('Normalized Firing Rate');
legend('Location','best');
box off;
% title(sprintf('PSTH (Array %s, Ch %d) — multi-target (t=0 Reach onset)', array_names(array), channel));


% ================== Raster ==================
ax2 = nexttile; hold on;

trial_counter = 0;

for tIdx = 1:numel(targets)
    tgt = targets(tIdx);
    this_color = colors_target(tgt,:);
    Trials_list = Trials_list_all{tIdx};

    for tr = 1:numel(Trials_list)
        trial_counter = trial_counter + 1;
        trial_spikes = Trials_list{tr};

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
                plot([st st], [trial_counter-0.2 trial_counter+0.2], ...
                     'Color', this_color);
            end
        end
    end

    % linea di separazione tra target diversi (opzionale)
    yline(trial_counter + 0.5, ':', 'Color', [0.5 0.5 0.5], ...
          'HandleVisibility','off');
end

xlabel('Time (s)');
ylabel('Trials');
box off;

% sincronizza gli assi del tempo
linkaxes([ax1 ax2],'x');
axis tight;
