% =========================================================
% DESCRIZIONE:
%  Script per eseguire una decodifica temporale del Target_ID
%  tramite SVM, utilizzando il firing rate medio estratto da
%  finestre temporali scorrevoli su dati neurali.
%
%  Funzionamento generale:
%    - viene caricato il dataset
%    - si definisce una finestra scorrevole di ampiezza w_length
%      con overlap pari al 50%
%
%  Per ciascuna finestra scorrevole:
%    - per ogni set e per ogni trial vengono estratti i dati
%      da entrambi gli array e concatenati lungo la dimensione
%      dei canali
%    - all’interno della finestra (N_w bin) vengono calcolati
%      i firing rate medi per ciascun canale (Hz)
%    - si costruisce così la matrice delle feature X
%      (tutti i trial × tutti i canali)
%    - si costruisce il vettore delle etichette Y (Target_ID)
%    - si esegue una k-fold cross-validation (k = 5) tramite
%      la funzione esterna svm_cv per addestrare un SVM multiclass
%      e stimare l’accuratezza della finestra corrente
%    - vengono salvate:
%         * l’accuratezza complessiva
%         * l’accuratezza per ciascuna classe (diag normalizzata)
%
%  Visualizzazione finale:
%    - andamento temporale dell’accuratezza complessiva e delle
%      accuratezze per classe (curve smussate nel tempo)
%    - indicazione degli eventi/stati del task tramite linee
%      verticali e relative etichette
%    - indicazione del livello di chance (1 / n_classes)
%    - asse temporale basato sul centro di ciascuna finestra
%
%  Parametri principali modificabili:
%    - n_sets, n_arrays, n_trials
%    - bin_size, w_length, overlap
%    - funzione svm_cv (tipo di SVM, k-fold, ecc.)
% =========================================================


clearvars
% close all
clc 

n_sets = 6;
n_arrays = 2;
n_trials = 32;
bin_size = 0.02;

filename = '../00_Data_extraction/gaze_BCI02.mat';
load(filename)

%% Decoding over time with SVM
% bins per trial
N = size(data(1).Data(1).Interp(1).Trial, 1); 
rec_duration = N*bin_size; 

% finestra scorrevole
w_length = 0.6; 
overlap = 0.5*w_length;
N_w = round(w_length/bin_size);
N_o = round(overlap/bin_size);

% etichette Y
Y = []; 
for set = 1:n_sets
    Y = [Y; [data(set).Data(1).Interp.Target_ID]']; 
end  
classes = unique(Y,'stable');
n_classes = numel(classes);

% numero finestre
n_acc = floor((N - N_w)/(N_w-N_o)) + 1; 
acc_overall = zeros(n_acc,1);
acc_class   = zeros(n_acc, n_classes);

% loop sulle finestre
start_w = 1; 
end_w = start_w + N_w - 1; 

for w = 1:n_acc
    X = cell(n_trials*n_sets,1);
    j = 1; 
    for set = 1:n_sets
        for trial = 1:n_trials
            SVM_matrix = []; 
            for array = 1:n_arrays
                SVM_matrix = [SVM_matrix, data(set).Data(array).Interp(trial).Trial(start_w:end_w, :)]; 
            end 
            X{j} = mean(SVM_matrix./bin_size,1);
            j = j + 1; 
        end   
    end 
    X = cell2mat(X); 
    
    % k-fold SVM
    k_fold = 5; 
    [acc_overall(w), cm{w}] = svm_cv(X, Y, k_fold);
    
    % accuracy per classe (diag della cm normalizzata per riga) 
    cm_norm = cm{w} ./ max(sum(cm{w},2),1);
    acc_class(w,:) = diag(cm_norm)*100;

    % prossima finestra
    start_w = start_w + (N_w - N_o); 
    end_w = start_w + N_w - 1;  
end 

%% Figure
figure('Color', 'White')
events_time_tmp = []; 
for i = 1:length(data(1).Data(2).Interp(1).Task_states)
    events_time = [events_time_tmp; size(data(1).Data(2).Interp(1).Task_states{i,2},1)*bin_size];
    events_time_tmp = events_time; 
end 
increment_times = cumsum(events_time); 
labels = string(data(1).Data(2).Interp(1).Task_states(:,1));
labels = ["", "Target cue", "Go cue", ""];

colors = [
    0.839, 0.153, 0.157;  % rosso
    0.122, 0.467, 0.706;  % blu
    0.172, 0.627, 0.172;  % verde
    0.580, 0.404, 0.741;  % viola
    1.000, 0.498, 0.055;  % arancione
    0.737, 0.741, 0.133;  % giallo oliva
    0.549, 0.337, 0.294;  % marrone
    0.890, 0.466, 0.760;  % rosa
];
alpha = 0.5;

t = (((0:n_acc-1)*(N_w - N_o)) + N_w/2)*bin_size;
w_smooth = 5;
for c = 1:n_classes
    acc_smooth = smoothdata(acc_class(:, c), 'gaussian', w_smooth);
    plot(t, acc_smooth, 'LineWidth', 1.0, 'Color', [colors(c,:),  alpha], 'HandleVisibility','off'), hold on
end
acc_smooth_overall = smoothdata(acc_overall, 'gaussian', w_smooth);
plot(t, acc_smooth_overall*100, 'LineWidth', 1.5, 'Color', 'k', 'DisplayName','Overall'), hold on

if exist('increment_times','var') && ~isempty(increment_times)
    xline(increment_times, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
    labels = {"Target cue", "Go cue"};

    ylim([0 100]);
    ax = gca;

    y_pos = ax.YLim(2) - 5;
    for i = 1:2
        x_pos = increment_times(i) - 0.3;
        text(x_pos - 0.05, y_pos, labels{i}, ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', ...
            'Rotation', 90, ...
            'FontSize', 11);
    end
end



yline((1/n_classes)*100,'-', 'Chance', 'HandleVisibility','off'); 
% legend show; 
xlabel('Time (s)');
ylabel('Accuracy (%)');
xlim([0 rec_duration]);
box off; 

