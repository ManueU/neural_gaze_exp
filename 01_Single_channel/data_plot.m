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

close all
clc

n_sets = 6; 
n_arrays = 2;
ch_start = 1;
ch_end = 96;
target_des = [1, 2, 3, 4, 5, 6, 7, 8]; 
bin_size = 0.02;

filename = "../00_Data_extraction/free-gaze_BCI02.mat";
load(filename);

%% Single channel
events_time_tmp = []; 
for i = 1:length(data(1).Data(2).Interp(1).Task_states)
    events_time = [events_time_tmp; size(data(1).Data(2).Interp(1).Task_states{i,2},1)*bin_size];
    events_time_tmp = events_time; 
end 
increment_times = cumsum(events_time); 

labels = string(data(1).Data(2).Interp(1).Task_states(:,1));

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

for array = 1:n_arrays
    for channel = ch_start:ch_end
        flag = 0; 

        for target = 1:length(target_des)
            M_spikes = [];

            % M_spikes contains the spike count matrix of each trial
            % associated to the specific target id. We consider all the
            % sets.
            for set = 1:n_sets
                idx = find([data(set).Data(array).Interp.Target_ID] == target_des(target));                      
                for j = 1:length(idx)
                    M_spikes = [M_spikes, [data(set).Data(array).Interp(idx(j)).Trial(:,channel)]];   
                end
            end 
            M_spikes_mean = mean(M_spikes, 2); 
            M_spikes_std  = std(M_spikes_mean, 0, 2);
            n_trials      = size(M_spikes, 2);
            M_spikes_sem  = M_spikes_std / sqrt(n_trials);
        
            firing_rate = M_spikes_mean ./ bin_size;
            firing_std  = M_spikes_std  ./ bin_size;            

            data_zscored = (firing_rate-mean_baseline(channel, array))./std_baseline(channel, array); 
            z_std = firing_std / std_baseline(channel, array);

        
            w = 25; 
            data_zscored_s = smoothdata(data_zscored, 'gaussian', w); 
            fr_s   = smoothdata(firing_rate, 'gaussian', w);
            
            if(max(fr_s) > 2)
                if flag == 0
                    figure('Color','w'); 
                    hold on
                    if exist('increment_times','var') && ~isempty(increment_times)
                       xline(increment_times, 'k', 'HandleVisibility','off');
                    end

                    xlabel('Time (s)');
                    ylabel('z-scored firing rate');
                    ylim([-10 10])
                    title(sprintf('Array = %s; Channel = %d;', array_names(array), channel));
                    legend('Location','best');
                    flag = 1; 
                end 
                
                t = (1:numel(firing_rate)) * bin_size;
                plot(t, data_zscored_s, 'LineWidth', 1.5, 'Color', colors_target(target,:), 'DisplayName', sprintf('%d', target_des(target))), hold on
                
                % yline(2, '--k', '2 std', 'LineWidth', 1, 'DisplayName', '2 std', 'HandleVisibility','off');
                % yline(-2, '--k', '-2 std', 'LineWidth', 1, 'DisplayName', '-2 std', 'HandleVisibility','off');

            end 
        end 

        if flag == 1
            ax = gca;
            x_times = [0; increment_times(:)];
            x_text  = x_times(1:end-1) + diff(x_times)/2;
        
            y_text = (ax.YLim(2) - 0.2) * ones(size(x_text));
            text(x_text, y_text, labels, 'HorizontalAlignment', 'center');
        end
    end 
end 

