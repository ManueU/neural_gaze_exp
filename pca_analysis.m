clear 
% close all
clc

%% Data extraction 1
[file, path] = uigetfile('*.mat', 'Seleziona un file MAT');
if isequal(file,0)
    error('Nessun file selezionato.');
else
    disp('Loading data...')
    load(fullfile(path, file));
end

MaskedSpikeCounts = set01.SpikeCount(:,1:5:end);
Motor1SC = MaskedSpikeCounts(:,[1:64 97:128]);
Stim1SC = MaskedSpikeCounts(:,65:96);
Motor2SC = MaskedSpikeCounts(:,[129:192 225:256]);
Stim2SC = MaskedSpikeCounts(:,193:224);

target_coordinates = [set01.TaskStateMasks.target(2, :)'  set01.TaskStateMasks.target(3, :)'];
target_info = target_coordinates; 
target_info(:,3) = NaN;
refXY = [  0      0
           0      0.2
           0.14   0.14
           0.2    0
           0.14  -0.14
           0     -0.2
          -0.14  -0.14
          -0.2    0
          -0.14   0.14 ];
codes = (0:8).';
tol = 1e-6;
[tf, loc] = ismembertol(target_info(:,1:2), refXY, tol, 'ByRows', true);
target_info(tf,3) = codes(loc(tf));

task_state = set01.TaskStateMasks.state_name'; 
trials = set01.trial_num'; 

n_targets = 8; 
pca_matrix = []; 
for channel = 1:96
    col_pca_matrix = [];
    for c = 1:n_targets
        idx = find(task_state == "Reach" & target_info(:,3) == c);
        d = diff(idx);
        ind_stacchi = [0 find(d > 1).' length(idx)];
        blocchi = cell(1, length(ind_stacchi)-1);
        for k = 1:length(ind_stacchi)-1
            blocchi{k} = idx(ind_stacchi(k)+1 : ind_stacchi(k+1));
        end
    
        len_blocchi = cellfun(@length, blocchi);
        len_min = min(len_blocchi); 
        disp(len_min);

        % Creo la matrice
        num_blocchi = numel(blocchi);
        M = zeros(len_min, num_blocchi);
        
        for k = 1:num_blocchi
            indici = blocchi{k}(1:len_min);      
            M(:,k) = Motor2SC(indici, channel).';         
        end
        firing_rate = mean(M,2)./0.02;
        col_pca_matrix = [col_pca_matrix; firing_rate]; 
    end 
    pca_matrix(:,channel) = col_pca_matrix; 
end


%% Data extraction 2
load('motor.mat')
n_trials = [32,32,32,32]; 

% Firing rate 
n_sets = 4; 
n_targets = 8; 
ch_start = 1; 
ch_end = 96; 
bin_size = 0.02; 

idx_reach = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == "Reach"); 
idx_prereach = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == "Pres12"); 

pca_matrix_array = []; 
for array = 2:2
    pca_matrix = []; 
    for channel = ch_start:ch_end
        firing_rate = []; 
        for target = 1:n_targets
            M_spikes_tmp = [];
            for set = 1:n_sets
                idx = find([data(set).Data(array).Resampled.Target_ID] == target);                      
                for j = 1:length(idx)
                    % M_spikes = [M_spikes_tmp, [dataset(set).Data(array).Resampled(idx(j)).Trial(:,channel)]];  
                    M_spikes = [M_spikes_tmp, [data(set).Data(array).Resampled(idx(j)).Task_states{idx_reach, 2}(1:70,channel)]];
                    M_spikes_tmp = M_spikes;
                end
            end 
            M_spikes_mean = mean(M_spikes, 2);     
            tmp = M_spikes_mean ./ bin_size;
            firing_rate = [firing_rate; tmp];
        end 
        pca_matrix = [pca_matrix, firing_rate]; 
    end
    pca_matrix_array = [pca_matrix_array, pca_matrix]; 
end 

%% pca
[Xz, muZ, sigmaZ] = zscore(pca_matrix, 0, 1); 
[coeff, score, latent, tsquared, explained] = pca(Xz, 'Algorithm','svd');

%%
mu = mean(pca_matrix);
sigma = std(pca_matrix);

i = 1; 
for c = 1 : n_targets
    data_c = pca_matrix(i:i+153, :);
    centered = (data_c - mu) ./ sigma;
    projected{c} = centered * coeff;
    i = i + 154; 
end


%% Figures
% 1
figure;
hold on;
target_ids = 1:8; 
for c = 1 : n_targets
    traj = projected{c};
    t1 = smoothdata(traj(:,1), 'gaussian', 1);
    t2 = smoothdata(traj(:,2), 'gaussian', 1);
    t3 = smoothdata(traj(:,3), 'gaussian', 1);

    plot3(t1, t2, t3,'LineWidth', 2, ...
        'DisplayName', num2str(target_ids(c)));
end
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
legend('Location', 'bestoutside');
grid on;
view(45, 25);


%% 2 
nbin = 70;       % numero di bin temporali
ntarget = 8;      % numero di target

figure('Color','White'); hold on
colors = [
    228 26 28;    % rosso
     55 126 184;  % blu
     77 175 74;   % verde
    152 78 163;   % viola
    255 127 0;    % arancione
    255 255 51;   % giallo
    166 86 40;    % marrone
    247 129 191;  % rosa
] / 255;
legend_entries = cell(ntarget,1);
% % % 
w = 10; 
for t = 1:ntarget
    idx = (t-1)*nbin + (1:nbin);   % righe di Xz/score che appartengono al target t
    traj = score(idx,1:3);         % traiettoria nei primi 3 PC
    t1 = smoothdata(traj(:,1), 'gaussian', w);
    t2 = smoothdata(traj(:,2), 'gaussian', w);
    t3 = smoothdata(traj(:,3), 'gaussian', w);

    plot3(t1, t2, t3, ...
          'Color', colors(t,:), 'LineWidth', 2, 'MarkerFaceColor', colors(t,:));
    legend_entries{t} = sprintf('Target %d', t);
end

xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
grid on; axis equal;
legend(legend_entries, 'Location', 'bestoutside'); 
% title('PCA');