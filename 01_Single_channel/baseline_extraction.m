clearvars
close all
clc

%% Baseline based on 60s quiet recording
[file, path] = uigetfile('*.mat', 'Seleziona un file MAT');
if isequal(file,0)
    error('Nessun file selezionato.');
else
    disp('Loading data...')
    load(fullfile(path, file));
end

counts = cellfun(@numel, baselineStruct.spikes.motor);
bas_m_medial = counts(1:96)./60;
bas_m_lateral = counts(97:192)./60;

baseline_const = [bas_m_medial, bas_m_lateral]; 

bin_size = 0.02; 
time_bins = 0:bin_size:60;
n_channels = 192;
for ch = 1:n_channels
    s = baselineStruct.spikes.motor{ch};
    tmp = histcounts(s, time_bins); 
    bin_baseline(: ,ch) = tmp/bin_size; 
end 

raw_des = length(dataset(1).Data(1).Resampled(1).Trial);
raw_act = size(bin_baseline,1);
nBlocks = floor(raw_act / raw_des);

for ch = 1:n_channels
    bin = bin_baseline(:,ch); 
    bin = bin(1 : nBlocks * raw_des);
    baseline(ch).mean =  mean(reshape(bin, raw_des, nBlocks), 2);
    baseline(ch).std = std(mean(reshape(bin, raw_des, nBlocks))); 
    baseline(ch).sem = std(mean(reshape(bin, raw_des, nBlocks)))/sqrt(length(mean(reshape(bin, raw_des, nBlocks)))); 
end 
clearvars -except dataset baseline baseline_const

%% Baseline as suggested by John
filename = "../00_Data_extraction/free-gaze_BCI02.mat";
load(filename);

% Single channel
n_sets = 6; 
n_arrays = 2;
n_channels = 96;
n_trials = 32;  
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

for array = 1:n_arrays
    for channel = 1:n_channels
        flag = 0; 
        M_spikes = [];
        for set = 1:n_sets
            for trial = 1:n_trials
                 M_spikes = [M_spikes; [data(set).Data(array).Resampled(trial).Trial(:,channel)]];   
            end
         end 
        M_spikes_mean = mean(M_spikes); 
        M_spikes_std  = std(M_spikes);
        M_spikes_sem  = std(M_spikes)/sqrt(length(M_spikes));
        
        firing_rate = M_spikes_mean ./ bin_size;
        firing_std  = M_spikes_std  ./ bin_size;  
        firing_sem  = M_spikes_sem  ./ bin_size;  

        mean_baseline(channel, array) = firing_rate;
        std_baseline(channel, array) = firing_std;
    
    end 
end 
