clear all
close all
clc 

period_pre = 0.1; 
period_reach = 0.5;
n_trials = 32; % trials per set 
bin_size = 0.02; 
 
mat_files = { ...
    'motor_BCI02.mat' ... 
    'controlled_BCI02.mat'
};

%% PCA 
for d = 1:numel(mat_files) 
    disp(mat_files(d)); 
    ds_name = mat_files{d};
    load(ds_name);

    % Matrix construction
    if strcmp(ds_name, 'controlled_BCI02.mat')
        PRES = "Gaze"; 
    else 
        PRES = "Pres12"; 
    end 
    idx_pres = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == PRES); 
    idx_reach = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == "Reach"); 
    
    start_pres = size(data(1).Data(2).Resampled(1).Task_states{idx_pres,2}, 1) - (period_pre/bin_size); 
    end_reach = period_reach/bin_size; 
    
    pca_matrix_array = [];
    for array = 1:2
        pca_matrix = []; 
        for set = 1:6
            for trial = 1:n_trials
                tmp_pres = data(set).Data(array).Resampled(trial).Task_states{idx_pres, 2}(start_pres:end, :); 
                tmp_reach = data(set).Data(array).Resampled(trial).Task_states{idx_reach,2}(1:end_reach, :); 
                m_trial = [tmp_pres; tmp_reach]; 
                
                pca_matrix =  [pca_matrix; m_trial ./ bin_size]; 
            end
        end 
        pca_matrix_array = [pca_matrix_array, pca_matrix]; 
    end
    
    if strcmp(ds_name, 'motor_BCI02.mat')
        x_motor = pca_matrix_array; 
    else
        x_free = pca_matrix_array; 
    end

end 
X = [x_free; x_motor];
X = zscore(X, 0, 1); 

[coeff, score, latent, ~, explained] = pca(X);

%% Figure
w = 10; 
traj = score(1:size(score,1)/2,1:3);        
t1 = smoothdata(traj(:,1), 'gaussian', w);
t2 = smoothdata(traj(:,2), 'gaussian', w);
t3 = smoothdata(traj(:,3), 'gaussian', w);

plot3(t1, t2, t3, ...
      'Color', 'b', 'LineWidth', 2);
hold on; 


traj = score(size(score,1)/2 + 1:end,1:3);        
t1 = smoothdata(traj(:,1), 'gaussian', w);
t2 = smoothdata(traj(:,2), 'gaussian', w);
t3 = smoothdata(traj(:,3), 'gaussian', w);

plot3(t1, t2, t3, ...
      'Color', 'r', 'LineWidth', 2);

xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
grid on;
