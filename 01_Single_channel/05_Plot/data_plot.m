% =========================================================
% DESCRIZIONE:
% Analizza l’attività di spike di più canali e array corticali per:
% - calcolare il firing rate medio per canale e target di movimento;
% - normalizzare il firing rate (z-score) rispetto a una baseline;
% - identificare i canali “attivi” (firing rate smussato > soglia);
% - visualizzare, per ciascun canale attivo, le traiettorie z-scorate
%   nel tempo per tutti i target, con colori diversi;
% - segnare sul grafico i cambi di stato del task e le relative etichette.
%
% I dati vengono estratti dalla struttura:
%   data(set).Data(array).Resampled(trial)
% utilizzando i campi Target_ID, Trial e Task_states.
%
% Prima dell’esecuzione, modificare i parametri principali:
%   filename      : percorso/nome del file .mat da caricare
%   n_sets        : numero di set/sedute da analizzare
%   n_arrays      : numero di array (es. “medial”, “lateral”)
%   target_des    : lista degli ID di target da considerare
%   ch_start/ch_end : intervallo di canali da analizzare
%   bin_size      : ampiezza del bin temporale (s)
%   mean_baseline, std_baseline : parametri di baseline per lo z-score
%   w             : ampiezza finestra di smoothing
%   soglia attività: valore minimo di firing rate smussato (es. 2 spk/s)
% =========================================================

% clearvars
% close all
clc

sets_plot = [1, 2, 3, 4, 5, 6];
n_sets = numel(sets_plot); 
n_arrays = 2;
n_trials = 16; 
n_channels = 96;
ch_start = 33;
ch_end = 33;
target_des = [1,2,3,4,5,6,7,8]; 
bin_size = 0.02;

filename = { 
    '../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat'};
    % , ...
    % '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat', ...
    % '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat'};


%% Single channel
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

for d = 1:numel(filename)
    load(filename{d}); 

    TS = data(1).Data(1).Interp(1).Task_states;   
    state_names = string(TS(:,1));
    state_dur_s = cellfun(@(x) size(x,1)*bin_size, TS(:,2));
    
    state_onset_s = [0; cumsum(state_dur_s(1:end-1))];  % onset di ogni stato
    state_end_s   = cumsum(state_dur_s);                % fine di ogni stato
    get_onset = @(name) state_onset_s( find(strcmpi(state_names, name), 1, 'first') );
    
    thisfile = string(filename{d}); 
    lines_t = [];
    lines_lab = strings(0,1);
    
    if contains(thisfile,"free-gaze","IgnoreCase",true)
        task_name = "Free-gaze";
        lines_t   = [ get_onset("pres12"); get_onset("reach") ];
        lines_lab = ["Target cue"; "Go cue"];
    
    elseif contains(thisfile,"motor","IgnoreCase",true)
        task_name = "Gaze-on-center";
        lines_t   = [ get_onset("pres12"); get_onset("reach") ];
        lines_lab = ["Target cue"; "Go cue"];
    
    elseif contains(thisfile,"controlled","IgnoreCase",true)
        task_name = "Gaze-on-target";
        lines_t   = [ get_onset("pres12"); get_onset("gaze"); get_onset("reach") ];
        lines_lab = ["Target cue"; "Go cue - gaze"; "Go cue"];
    
    elseif contains(thisfile,"gaze","IgnoreCase",true)
        task_name = "Gaze-only";
        lines_t   = [ get_onset("pres12"); get_onset("gaze") ];
        lines_lab = ["Target cue"; "Go cue - gaze"];
    end


    for array = 1:n_arrays
        for channel = ch_start:ch_end
            ch_global = (array-1)*n_channels + channel;
            disp(ch_global);
            
            flag = 0; 
            for target = 1:length(target_des)
                firing_rate = []; 
    
                for set_ = 1:n_sets
                    set = sets_plot(set_);
                    idx = find([data(set).Data(array).Interp.Target_ID] == target_des(target) & [data(set).Data(array).Interp.Excluded] == 0);                      
                    for j = 1:length(idx)
                        M_spikes = data(set).Data(array).Interp(idx(j)).Trial(:,channel);  
                        firing_rate = [firing_rate, M_spikes ./ bin_size]; 
                    end
                end         
                % data_zscored = (mean(firing_rate, 2)-mean_baseline_common(ch_global)) ./ std_baseline_common(ch_global); 

                w = 15; 
                data_zscored_s = smoothdata(mean(firing_rate, 2), 'gaussian', w); 
                if flag == 0
                    figure('Color','w'); 
                    hold on
                    ax = gca;
                    for k = 1:numel(lines_t)
                        xline(lines_t(k), '--k', 'HandleVisibility','off');
                    end
                    % yline(-2)
                    % yline(2)


                    xlabel('Time (s)');
                    ylabel('Firing rate');
                    title({['\bf' char(task_name)], ['\rmChannel ' num2str(ch_global)]}, 'Interpreter','tex');
                    % legend('Location','best');
                    flag = 1; 
                end 
                t = (1:numel(data_zscored_s)) * bin_size;
                plot(t, data_zscored_s, 'LineWidth', 1.5, 'Color', colors_target(target,:), 'DisplayName', sprintf('%d', target_des(target))), hold on

            end 
            axis tight;
            dx = 0.05 * range(ax.XLim);          
            yl = ax.YLim;
            y_txt = yl(2) - 0.05*range(yl);      
            
            text(lines_t - dx, y_txt * ones(size(lines_t)), lines_lab, 'HorizontalAlignment','right', 'VerticalAlignment','top', 'Rotation',90, 'Clipping','on', 'HandleVisibility','off');
 
        end 
    end 
end

