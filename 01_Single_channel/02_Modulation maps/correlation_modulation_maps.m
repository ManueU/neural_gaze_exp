clearvars -except mean_baseline_common std_baseline_common
close all
clc

%% mappa somatotopica inspired by Natalya's paper
filename_all = {'../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat', ...
                '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat', ...
                '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat'};


corr_data = struct(); 
corr_data.n_condition = 3;
corr_data.n_targets = n_targets;
corr_data.array = 2;

labels = {"Free-gaze", "Gaze-on-center", "Gaze-on-target"};
for d = 1:numel(filename_all)
    filename = filename_all{d}; 

    [~, baseName, ext] = fileparts(filename);
    ds_name = [baseName ext];

    if strcmp(ds_name, 'controlled_BCI02_exclUpdated.mat')
        PRE = "Gaze";
        POST = "Reach";
        file = 'responsive_channels_gaze_on_target.mat'; 
    else
        PRE = "Pres12";
        POST = "Reach";
        if strcmp(ds_name, 'motor_BCI02_exclUpdated.mat')
            file = 'responsive_channels_gaze_on_center.mat';
        else
            file = 'responsive_channels_free_gaze.mat';
        end 
    end

    run('modulation_mask_analysis.m');
    run('modulation_maps_analysis.m');

    motor_map = motor_electrodes{corr_data.array};
    valid_outline = ~isnan(motor_map);
    for target = 1:n_targets
        masked_index = modulation_matrix{target, corr_data.array} .* modulation_mask{1,corr_data.array};
        smoothed_modulation_matrix{target} = nangauss_smooth(masked_index, 4, 0.5);
    end
    
    correlation_maps{d} = gradientTargetPreferenceCircular(smoothed_modulation_matrix, corr_data, labels, d, valid_outline);
end

%% Correlation between maps
correlationBetweenTargetPreferenceMaps(correlation_maps, corr_data, labels)
