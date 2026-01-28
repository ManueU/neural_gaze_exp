%% Two-way ANOVA
clearvars -except mean_baseline std_baseline
% close all
clc

%% PARAMETRI
sets = [2,4,5,6]; 
n_sets = numel(sets); 
n_arrays = 2;      
n_channels = 96;    
n_targets = 8; 
bin_size = 0.02; 
period_pre = 0.1;   
period_post = 0.5;   

filename = { ...
    '../00_Data_extraction/free-gaze_BCI02.mat', ...
    '../00_Data_extraction/motor_BCI02.mat',...
    '../00_Data_extraction/controlled_BCI02.mat'
    };

% Etichette per il fattore gaze (stesso ordine di filename)
gaze_label = {'free', 'center', 'target'};
all_target = [];   % target ID per trial
all_gaze = [];     % condizione di gaze per trial (1,2,3)

%% LOOP SU CONDIZIONI (GAZE)
condition = [];
for d = 1:numel(filename)
    fprintf('\nDataset: %s  (gaze = %s)\n', filename{d}, gaze_label{d}); 
    load(filename{d}); 
    
    [~, baseName, ext] = fileparts(filename{d});
    ds_name = [baseName ext];

    if strcmp(ds_name, 'controlled_BCI02.mat')
        PRE = "Gaze";
        POST = "Reach";
    elseif strcmp(ds_name, 'gaze_BCI02.mat')
        PRE = "Pres12";
        POST = "Gaze";
    else
        PRE = "Pres12";
        POST = "Reach";
    end
        
    idx_pres  = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == PRE); 
    idx_reach = find(string(data(1).Data(1).Interp(1).Task_states(:,1)) == POST); 
    if isempty(idx_pres) || isempty(idx_reach)
        error('Stati PRE/POST non trovati: controlla PRE/POST');
    end
    
    all_FR = [];
    pre_bins = max(1, round(period_pre/bin_size));
    post_bins = max(1, round(period_post/bin_size));

    fr_arrays = [];
    for array = 1:n_arrays
        ch_global = (array-1)*n_channels + n_channels;

        fr_trial = []; 
        trial_idx = 1;
        for set_ = 1:n_sets
            set = sets(set_);
            % Numero di trial in questo set/array
            n_trials_set = numel(data(set).Data(array).Interp);

            for tr = 1:n_trials_set
                thisTrial = data(set).Data(array).Interp(tr);
                targ = thisTrial.Target_ID;

                pres_state  = thisTrial.Task_states{idx_pres,  2}; 
                reach_state = thisTrial.Task_states{idx_reach, 2}; 
                
                tmp_pre  = pres_state(end-pre_bins+1:end, :); 
                tmp_post = reach_state(1:post_bins, :); 
                
                window_spikes = [tmp_pre; tmp_post]./bin_size;   
                norm_fr = (window_spikes - mean_baseline{1,d}(trial_idx,ch_global-n_channels + 1:ch_global)); 
                fr_trial = [fr_trial; mean(norm_fr,1)];
                if array == 1
                    all_target = [all_target; targ];
                    all_gaze   = [all_gaze; d];        % 1=free,2=center,3=target
                end 
                trial_idx = trial_idx + 1;

            end
        end 
        fr_arrays = [fr_arrays, fr_trial];  
     end
     condition = [condition; fr_arrays];
end
fprintf('\nTotale trial raccolti: %d\n', size(condition,1));

%% ANOVA 2-way (target × gaze) per ogni canale
n_tot_channels = n_arrays * n_channels;

p_target = nan(n_tot_channels,1);
p_gaze = nan(n_tot_channels,1);
p_interact = nan(n_tot_channels,1);

for ch = 1:n_tot_channels
    y = condition(:, ch);   
    
    [p, ~, stats] = anovan(y, {categorical(all_target), categorical(all_gaze)}, 'model', 'interaction', 'varnames', {'Target','Gaze'}, 'display','off');
    p_target(ch) = p(1);
    p_gaze(ch) = p(2);
    p_interact(ch) = p(3);

end

%% Correzione multipla FDR Benjamini-Hochberg
[h_target, ~, ~, adjp_target] = fdr_bh(p_target);
[h_gaze,   ~, ~, adjp_gaze]   = fdr_bh(p_gaze);
[h_inter,  ~, ~, adjp_inter]  = fdr_bh(p_interact);

alpha = 0.01; 
is_target = adjp_target < alpha;
is_gaze = adjp_gaze < alpha;
is_interact = adjp_inter < alpha;

fprintf('\n--- Risultati ANOVA ---\n');
fprintf('Canali mod. da TARGET:   %d / %d\n', sum(is_target), n_tot_channels);
fprintf('Canali mod. da GAZE:     %d / %d\n', sum(is_gaze), n_tot_channels);
fprintf('Canali con INTERAZIONE:  %d / %d\n', sum(is_interact), n_tot_channels);

%% Classificazione canali
% 0 = non significativo
% 1 = target-only
% 2 = gaze-only
% 3 = both (mixed)
% 4 = interaction-only (senza main effects)

chan_class = zeros(n_tot_channels,1);

for ch = 1:n_tot_channels
    if is_target(ch) && ~is_gaze(ch) && ~is_interact(ch)
        chan_class(ch) = 1;
    elseif ~is_target(ch) && is_gaze(ch) && ~is_interact(ch)
        chan_class(ch) = 2;
    elseif is_target(ch) && is_gaze(ch)
        chan_class(ch) = 3;
    elseif ~is_target(ch) && ~is_gaze(ch) && is_interact(ch)
        chan_class(ch) = 4;
    else
        chan_class(ch) = 0;
    end
end

fprintf('\nClassi canali:\n');
fprintf('Target-only:   %d\n', sum(chan_class==1));
fprintf('Gaze-only:     %d\n', sum(chan_class==2));
fprintf('Mixed (both):  %d\n', sum(chan_class==3));
fprintf('Interaction:   %d\n', sum(chan_class==4));
fprintf('Non sig.:      %d\n', sum(chan_class==0));

%% Figure (1)

categories = {'Target','Gaze','Non-sig'};
values = [
    sum(chan_class==1), ...
    sum(chan_class==2), ...
    sum(chan_class==0)
];

percent = 100 * values / n_tot_channels;

figure('Color', 'w');
bar(percent, 'FaceColor',[0.2 0.4 0.8]);
clear set
set(gca,'XTickLabel',categories,'XTickLabelRotation',15,'FontSize',12);
ylabel('Percentage of channels (%)');
ylim([0 100]);
grid on;
box off; 


%% Mappa canali basata sull'ANOVA (chan\_class)
load('ChannelMap_BCI02.mat'); 
motor_medial      = ChannelMap.ChannelNumbers{1,1};
motor_lateral     = ChannelMap.ChannelNumbers{1,3};
motor_electrodes  = {motor_medial, motor_lateral};

class_mask = cell(1, n_arrays);

for array = 1:n_arrays
    M = motor_electrodes{array};   

    [n_rows, n_cols] = size(M);
    mask_array = nan(n_rows, n_cols);

    for i = 1:n_rows
        for j = 1:n_cols
            ch = M(i,j);    % numero di canale "globale" (1..192)

            if isnan(ch)
                continue;
            end

            % Usa direttamente la classe dell'ANOVA:
            % 0 = non sig, 1 = target, 2 = gaze, 3 = both, 4 = interaction
            mask_array(i,j) = chan_class(ch);
        end
    end

    class_mask{array} = mask_array;
end

%% Figure: mappe array basate su ANOVA

for array = 1:n_arrays
    
    M_elec = motor_electrodes{array};
    M_class = class_mask{array};

    [n_rows, n_cols] = size(M_elec);

    figure('Color','w'); 
    hold on;
    axis equal;
    axis ij;  
    xlim([0.5, n_cols + 0.5]);
    ylim([0.5, n_rows + 0.5]);
    clear set
    set(gca, 'XTick', 1:n_cols, 'YTick', 1:n_rows, 'FontSize', 10);
    box on;

    % --- Palette (5 colori chiari) ---
    c_gray    = [0.85 0.85 0.85];   % grigio chiaro
    c_blue    = [0.30 0.55 0.85];   % blu medio
    c_teal    = [0.20 0.70 0.70];   % teal / turchese
    c_orange  = [0.95 0.65 0.30];   % arancio caldo
    c_purple  = [0.65 0.45 0.85];   % viola medio
    
    % Mappa colori per classi 0..4
    classColor = cell(5,1);
    classColor{1} = c_gray;     % 0 non significativo
    classColor{2} = c_blue;     % 1 target-only
    classColor{3} = c_teal;     % 2 gaze-only
    classColor{4} = c_orange;   % 3 mixed
    classColor{5} = c_purple;   % 4 interaction-only


    title(sprintf('Array %d', array));
    for i = 1:n_rows
        for j = 1:n_cols

            ch = M_elec(i,j);

            if isnan(ch)
                continue;
            end

            this_class = M_class(i,j);
            idx = this_class + 1;          % this_class in [0..4]
            faceColor = classColor{idx};


            rectangle('Position',[j-0.5, i-0.5, 1, 1], 'FaceColor', faceColor, 'EdgeColor', 'None');
            text(j, i, num2str(ch), 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize',8);
        end
    end

    hold on;
    labels = { ...
    'No modulated', ...
    'Modulated by target', ...
    'Modulated by gaze', ...
    'Target + Gaze (mixed)', ...
    'Target \times Gaze (interaction)' ...
    };
    hLeg = gobjects(5,1);
    for k = 1:5
        hLeg(k) = patch(NaN, NaN, classColor{k}, 'EdgeColor', 'none');
    end
    legend(hLeg, labels, 'Location', 'northeastoutside', 'Box', 'off');
end

