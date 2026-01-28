clearvars -except mean_baseline std_baseline
close all
clc

%% mappa somatotopica inspired by Natalya's paper
filename_all = {'../00_Data_extraction/free-gaze_BCI02.mat', ...
            '../00_Data_extraction/motor_BCI02.mat', ...
            '../00_Data_extraction/controlled_BCI02.mat'};

filename_excel_all = {'visual_analysis_FG.xlsx', ...
                  'visual_analysis_GC.xlsx', ...
                  'visual_analysis_GT.xlsx'};

corr_data.n_condition = 3;
corr_data.n_targets = 8;
corr_data.array = 2;
labels = {"Free-gaze", "Gaze-on-center", "Gaze-on-target"};
for d = 1:numel(filename_all)
    filename = filename_all{d}; 
    filename_excel = filename_excel_all{d};

    [~, baseName, ext] = fileparts(filename);
    ds_name = [baseName ext];

    if strcmp(ds_name, 'controlled_BCI02.mat')
        PRE = "Gaze";
        POST = "Reach";
    else
        PRE = "Pres12";
        POST = "Reach";
    end

    run('modulation_mask_analysis.m');
    run('modulation_maps_analysis.m');

    for target = 1:n_targets
        masked_index = modulation_matrix{target, corr_data.array} .* modulation_mask{1,corr_data.array};
        smoothed_modulation_matrix{target} = nangauss_smooth(masked_index, 4, 0.5);
    end

    correlation_maps{d} = gradientTargetPreference(smoothed_modulation_matrix, corr_data, labels, d);
end

%% Correlation between maps
correlationBetweenTargetPreferenceMaps(correlation_maps, corr_data, labels)
