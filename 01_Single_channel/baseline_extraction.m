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

mean_baseline = nan(n_channels, n_arrays);
std_baseline  = nan(n_channels, n_arrays);
for array = 1:n_arrays
    for channel = 1:n_channels
        trial_baselines = [];
        for set = 1:n_sets
            for trial = 1:n_trials
                 trial_spikes = data(set).Data(array).Interp(trial).Trial(:,channel);  
                 trial_baseline_counts = mean(trial_spikes); 
                 trial_baselines(end+1) = trial_baseline_counts;            
            end
        end 
        baseline_mean_counts = mean(trial_baselines);
        baseline_std_counts  = std(trial_baselines);

        mean_baseline(channel, array) = baseline_mean_counts / bin_size; 
        std_baseline(channel, array)  = baseline_std_counts  / bin_size;
    end 
end 
