close all 
clc

[file, path] = uigetfile('*.mat', 'Seleziona un file MAT');
if isequal(file,0)
    error('Nessun file selezionato.');
else
    disp('Loading data...')
    load(fullfile(path, file));
end

%% Baseline
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
clearvars -except dataset baseline 


