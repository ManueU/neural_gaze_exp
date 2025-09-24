clear 
close all
clc

load('free-gaze.mat')
array = 2; 
n_trials = [40, 32, 32, 32]; 

%% Data extraction
% Target IDs
n_sets = 4; 
tmp = []; 
for set = 1:n_sets
    Y = [tmp; [dataset(set).Data(array).Resampled.Target_ID]']; 
    tmp = Y;
end 

% Firing rate 
n_target = 8; 
ch_start = 1; 
ch_end = 96; 
bin_size = 0.02; 
idx_reach = find(string(dataset(1).Data(1).Resampled(1).Task_states(:,1)) == "Reach"); 

pca_matrix_array = []; 
for array = 1:2
    pca_matrix = []; 
    for channel = ch_start:ch_end
        firing_rate = []; 
        for target = 1:n_target
            M_spikes_tmp = [];
            for set = 1:n_sets
                idx = find([dataset(set).Data(array).Resampled.Target_ID] == target);                      
                for j = 1:length(idx)
                    % M_spikes = [M_spikes_tmp, [dataset(set).Data(array).Resampled(idx(j)).Trial(:,channel)]];  
                    M_spikes = [M_spikes_tmp, [dataset(set).Data(array).Resampled(idx(j)).Task_states{idx_reach, 2}(:,channel)]];
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


% Data preprocessing
Xz = zscore(pca_matrix, 0, 1); 

% pca
[coeff, score, latent, tsquared, explained, mu] = pca(Xz, 'Algorithm','svd');

%%
mu = mean(pca_matrix);
sigma = std(pca_matrix);

i = 1; 
for c = 1 : n_target
    data_c = pca_matrix(i:i+154, :);
    centered = (data_c - mu) ./ sigma;
    projected{c} = centered * coeff;
    i = i + 155; 
end


%% Figures
figure;
hold on;
target_ids = 1:8; 
for c = [1 : n_target] % 10]% 5 11 : 15] % 1 : num_cond
    traj = projected{c};
    t1 = smoothdata(traj(:,1), 'gaussian', 10);
    t2 = smoothdata(traj(:,2), 'gaussian', 10);
    t3 = smoothdata(traj(:,3), 'gaussian', 10);

    plot3(t1, t2, t3,'LineWidth', 2, ...
        'DisplayName', num2str(target_ids(c)));
end
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
legend('Location', 'bestoutside');
grid on;
view(45, 25);