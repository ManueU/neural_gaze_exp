clear all
close all
clc

[file, path] = uigetfile('*.mat', 'Seleziona un file MAT');
if isequal(file,0)
    error('Nessun file selezionato.');
else
    disp('Loading data...')
    load(fullfile(path, file));
end

%%
clearvars -except allData
close all
clc

%% Lateral array
% Spike count

SpikeCount_matrix = set08.SpikeCount(:, 1:5:end);
% SpikeCount_matrix = allData.set08.data.SpikeCount(:, 1:5:end);
idx.medial_sens   = 65:96;                   % Sensory medial (32ch)
idx.lateral_sens  = 193:224;                 % Sensory lateral (32ch)
idx.medial_motor  = [1:64,   97:128];        % Motor medial (96ch)
idx.lateral_motor = [129:192, 225:256];      % Motor lateral (96ch)

arrays.motor_m = SpikeCount_matrix(:, idx.medial_motor);
arrays.motor_l = SpikeCount_matrix(:, idx.lateral_motor);
arrays.sens_m  = SpikeCount_matrix(:, idx.medial_sens);
arrays.sens_l  = SpikeCount_matrix(:, idx.lateral_sens);

% Target info
target_coordinates = set08.TaskStateMasks.target(2:3, :);
% target_coordinates = allData.set08.data.TaskStateMasks.target(2:3, :);

% Label trial
bins_per_trial = countcats(categorical(set08.trial_num));
% bins_per_trial = countcats(categorical(allData.set08.data.trial_num));
idx_end_trial = cumsum(bins_per_trial);

trials.motor_l = cell(1, length(idx_end_trial));
trials.sens_l = cell(1, length(idx_end_trial));
for i = 1:length(idx_end_trial)
    if i == 1
        start_idx = 1;
    else
        start_idx = idx_end_trial(i-1) + 1;
    end

    end_idx = idx_end_trial(i);
    trials.motor_l{i} = arrays.motor_l(start_idx:end_idx, :);
    trials.sens_l{i} = arrays.sens_l(start_idx:end_idx, :);
    target_per_bins{i} = target_coordinates(:,start_idx:end_idx);
end

nRows = cellfun(@(x) size(x,1), trials.motor_l);
minLen = min(nRows);
for i = 1:numel(trials.motor_l)
    trials.motor_l{i} = trials.motor_l{i}(1:minLen, :);
end

nRows = cellfun(@(x) size(x,1), trials.sens_l);
minLen = min(nRows);
for i = 1:numel(trials.sens_l)
    trials.sens_l{i} = trials.sens_l{i}(1:minLen, :);
end

% Label target
% XY_des = [ 0      0.2
%            0.2    0
%            0     -0.2
%           -0.2    0];

XY_des = [  0      0.2
           0.14   0.14
           0.2    0
           0.14  -0.14
           0     -0.2
          -0.14  -0.14
          -0.2    0
          -0.14   0.14 
          ];

n_trials = 32; 
target_per_trial = nan(1,n_trials); 

for i = 1:n_trials
    [tf, which] = ismember(target_per_bins{1,i}(1:2,:).', XY_des, 'rows');  
    target_per_trial(i) = max(which);
end

%% Plot sensory and motor - PSTH
n_targets = 8; 
n_ch_motor = 96;
n_ch_sens = 32;
ch_start = 1;
bin_size = 0.02;
w = 15;

% Motor
for channel = ch_start:n_ch_motor
    figure()
    for target = 1:n_targets 
        ch_fr = [];
        idx = find(target_per_trial == target);
        for j = 1:length(idx)
            ch_fr = [ch_fr, trials.motor_l{1,idx(j)}(:,channel)/bin_size];
        end
        ch_fr_s = smoothdata(mean(ch_fr,2), 'gaussian', w);
        plot(ch_fr_s)
        leg(target) = "Target " + target;
        title("Motor lateral")
        hold on 
    end
    legend(leg, 'Location', 'best')
    hold off
end

% Sensory
for channel = 1:n_ch_sens
    figure(channel)
    for target = 1:n_targets 
        ch_fr = [];
        idx = find(target_per_trial == target);
        for j = 1:length(idx)
            ch_fr = [ch_fr, trials.sens_l{1,idx(j)}(:,channel)/bin_size];
        end
        ch_fr_s = smoothdata(mean(ch_fr,2), 'gaussian', w);
        plot(ch_fr_s)
        title("Sensory lateral")
        hold on 
    end
    hold off
end

%% PSTH + Raster - Motor 
ch_start = 1;
ch_end = 30;

nBins = size(trials.motor_l{1,1}, 1);
t = ((1:nBins) - 0.5) * bin_size;

for channel = ch_start:ch_end
    for target = 1:1

        idx = find(target_per_trial == target);
        nTr = numel(idx);

        ch_counts = zeros(nBins, nTr);
        for j = 1:nTr
            ch_counts(:, j) = trials.motor_l{1, idx(j)}(:, channel);
        end

        fr_mean = mean(ch_counts, 2) / bin_size;        
        fr_s    = smoothdata(fr_mean, 'gaussian', w);

        figure('Name', sprintf('Motor L - Ch %d - Target %d', channel, target)); clf
        tiledlayout(2,1, "TileSpacing","compact", "Padding","compact");

        ax1 = nexttile;
        plot(t, fr_s, 'LineWidth', 1.5);
        title(sprintf('Motor lateral | Channel %d | Target %d', channel, target))
        ylabel('FR (Hz)')

        ax2 = nexttile;
        hold on
        
        for j = 1:nTr
            counts = ch_counts(:, j);
        
            for b = 1:nBins
                k = counts(b);
                if k <= 0, continue; end
        
                % limiti del bin (in secondi)
                t0 = (b-1) * bin_size;
                t1 = b * bin_size;
        
                % k spike -> k posizioni dentro il bin (evito bordi con k+1)
                spk_t = t0 + (1:k) * (bin_size/(k+1));
        
                % disegno i tick
                for s = 1:numel(spk_t)
                    plot([spk_t(s) spk_t(s)], [j-0.4 j+0.4], 'k', 'LineWidth', 0.5);
                end
            end
        end
        
        ylim([0 nTr+1])
        xlim([0 nBins*bin_size])
        ylabel('Trial')
        xlabel('Time (s)')
        box on
        hold off
    end
end
