% =========================================================
% DESCRIZIONE:
%  Script per eseguire una decodifica temporale del Target_ID
%  tramite SVM, utilizzando il firing rate medio estratto da
%  finestre temporali scorrevoli su dati neurali.
%
%    - decoding con 4-fold repeated cross-validation
%    - loop automatico su 3 condizioni:
%         * free-gaze
%         * motor
%         * controlled
%    - asse temporale allineato al Go cue / Go cue cursor
%    - correzione del posizionamento delle etichette evento:
%      gli eventi pre-0 vengono etichettati a destra della linea,
%      quelli post-0 a sinistra, così il testo non esce dal pannello
% =========================================================

clearvars
close all
clc

n_sets   = 6;
n_arrays = 2;
n_trials = 16;
bin_size = 0.02;

filenames = { ...
    '../00_Data_extraction/BCI02_Session_0924/free-gaze_BCI02_exclUpdated.mat'
    % , ...
    % '../00_Data_extraction/BCI02_Session_0924/motor_BCI02_exclUpdated.mat', ...
    % '../00_Data_extraction/BCI02_Session_0924/controlled_BCI02_exclUpdated.mat'
    };

cond_names = {'Free-gaze', 'Gaze-on-center', 'Gaze-on-target'};

all_results = cell(numel(filenames),1);

for d = 1:numel(filenames)

    fprintf('\n=================================================\n');
    fprintf('Running condition: %s\n', cond_names{d});
    fprintf('File: %s\n', filenames{d});
    fprintf('=================================================\n');

    load(filenames{d})

    %% Decoding over time with SVM
    N = size(data(1).Data(1).Interp(1).Trial, 1);
    rec_duration = N * bin_size;

    % finestra scorrevole
    w_length = 0.6;
    overlap  = 0.5 * w_length;
    N_w = round(w_length / bin_size);
    N_o = round(overlap  / bin_size);

    % etichette Y
    Y = [];
    for set = 1:n_sets
        idx = [data(set).Data(1).Interp.Excluded] == 0;
        Y = [Y; [data(set).Data(1).Interp(idx).Target_ID]'];
    end

    classes   = unique(Y,'stable');
    n_classes = numel(classes);

    % numero finestre
    n_acc = floor((N - N_w) / (N_w - N_o)) + 1;

    acc_overall    = zeros(n_acc,1);
    balacc_overall = zeros(n_acc,1);
    recall_class   = zeros(n_acc, n_classes);

    acc_std    = zeros(n_acc,1);
    balacc_std = zeros(n_acc,1);

    start_w = 1;
    end_w   = start_w + N_w - 1;

    n_valid = sum(arrayfun(@(s) sum([data(s).Data(1).Interp.Excluded] == 0), 1:n_sets));

    % parametri CV
    Ycat     = categorical(Y);
    counts   = countcats(Ycat);
    minCount = min(counts);

    k_fold    = min(6, minCount);
    n_repeats = 20;
    seed      = 0;

    if k_fold < 2
        error('Una o più classi hanno <2 campioni: impossibile fare CV.');
    end

    fprintf('Class counts: ');
    fprintf('%d ', counts);
    fprintf('\n');
    fprintf('Using repeated %d-fold CV with %d repeats\n', k_fold, n_repeats);

    cm = cell(n_acc,1);
    metrics = cell(n_acc,1);

    % loop sulle finestre
    for w = 1:n_acc

        X = cell(n_valid, 1);
        j = 1;

        for set = 1:n_sets
            for trial = 1:n_trials
                SVM_matrix = [];

                if data(set).Data(1).Interp(trial).Excluded == 0
                    for array = 1:n_arrays
                        SVM_matrix = [SVM_matrix, data(set).Data(array).Interp(trial).Trial(start_w:end_w, :)];
                    end

                    X{j} = mean(SVM_matrix ./ bin_size, 1);
                    j = j + 1;
                end
            end
        end

        X = cell2mat(X);

        % repeated 4-fold SVM
        [results, cm_pooled] = svm_repeated_cv(X, Y, k_fold, n_repeats, seed);

        cm{w} = cm_pooled;
        metrics{w} = results;

        acc_overall(w)    = mean(results.acc);
        balacc_overall(w) = mean(results.balAcc);

        acc_std(w)    = std(results.acc);
        balacc_std(w) = std(results.balAcc);

        recall_class(w,:) = mean(results.perClassRecall, 1) * 100;

        start_w = start_w + (N_w - N_o);
        end_w   = start_w + N_w - 1;
        
    end

    %% Figure
    figure('Color','w');
    hold on;

    % tempi degli eventi
    events_time_tmp = [];
    for i = 1:length(data(1).Data(2).Interp(1).Task_states)
        events_time = [events_time_tmp; size(data(1).Data(2).Interp(1).Task_states{i,2},1) * bin_size];
        events_time_tmp = events_time;
    end
    increment_times = cumsum(events_time);

    % centro delle finestre
    t = (((0:n_acc-1) * (N_w - N_o)) + N_w/2) * bin_size;

    % allineamento temporale al Go cue / Go cue cursor
    switch d
        case 1   % Free-gaze
            align_time = increment_times(2);   % Go cue
        case 2   % Motor
            align_time = increment_times(3);   % Go cue
        case 3   % Controlled
            align_time = increment_times(4);   % Go cue cursor
    end

    t = t - align_time;
    increment_times = increment_times - align_time;

    % ===== VERIFICA finestra equivalente a [-0.1, +0.5] =====
    [~, idx_cmp] = min(abs(t - 0.2));
    fprintf('\nCondition: %s\n', cond_names{d});
    fprintf('Closest window center to 0.2 s: %.3f s\n', t(idx_cmp));
    
    t_start = t(idx_cmp) - w_length/2;
    t_end   = t(idx_cmp) + w_length/2;
    fprintf('Window limits: [%.3f, %.3f] s\n', t_start, t_end);
    fprintf('Window index: %d\n', idx_cmp);

    % smoothing
    w_smooth   = 5;
    acc_smooth = smoothdata(acc_overall * 100,    'gaussian', w_smooth);
    bal_smooth = smoothdata(balacc_overall * 100, 'gaussian', w_smooth);

    recall_class_smooth = zeros(size(recall_class));
    for c = 1:n_classes
        recall_class_smooth(:,c) = smoothdata(recall_class(:,c), 'gaussian', w_smooth);
    end

    colors = [
        0.839, 0.153, 0.157;
        0.122, 0.467, 0.706;
        0.172, 0.627, 0.172;
        0.580, 0.404, 0.741;
        1.000, 0.498, 0.055;
        0.737, 0.741, 0.133;
        0.549, 0.337, 0.294;
        0.890, 0.466, 0.760;
    ];

    if n_classes > size(colors,1)
        colors = repmat(colors, ceil(n_classes / size(colors,1)), 1);
    end

    % recall per classe
    for c = 1:n_classes
        plot(t, recall_class_smooth(:,c), ...
            'LineWidth', 0.9, ...
            'Color', colors(c,:), ...
            'HandleVisibility','off');
    end

    % curve globali
    p1 = plot(t, acc_smooth, 'k',   'LineWidth', 1.5, 'DisplayName','Accuracy');
    p2 = plot(t, bal_smooth, 'k--', 'LineWidth', 1.5, 'DisplayName','Balanced accuracy');

    % chance level
    chance_bal = (1 / n_classes) * 100;
    yline(chance_bal, ':', 'Chance', 'HandleVisibility','off');

    % label events
    if ~isempty(increment_times)

        switch d
            case 1   % Free-gaze
                labels = {'Target cue', 'Go cue'};
                event_idx = 1:2;

            case 2   % Motor
                labels = {'', 'Target cue', 'Go cue'};
                event_idx = 2:3;

            case 3   % Controlled
                labels = {'', 'Target cue', 'Go cue - gaze', 'Go cue - cursor'};
                event_idx = 2:4;
        end

        % linee evento
        xline(increment_times(event_idx), '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');

        ylim([0 100]);
        ax = gca;
        y_pos = ax.YLim(2) - 10;


        for i = event_idx

            if isempty(labels{i}) || strcmp(labels{i}, "")
                continue
            end

            x_text = increment_times(i) + 0.05;
            h_align = 'right';

            text(x_text, y_pos, labels{i}, ...
                'HorizontalAlignment', h_align, ...
                'VerticalAlignment', 'top', ...
                'Rotation', 90, ...
                'FontSize', 11);
        end
    end

    ylim([0 100]);
    xlim([min(t)-w_length/2 max(t)]);
    yticks(0:20:100);
    xlabel('Time (s)');
    ylabel('Performance (%)');
    title(sprintf('%s' , cond_names{d}));
    box off;
    legend([p1 p2], 'Location', 'best');
    grid on

    % salvataggio risultati
    all_results{d} = struct( ...
        'condition', cond_names{d}, ...
        'filename', filenames{d}, ...
        't', t, ...
        'acc_overall', acc_overall, ...
        'balacc_overall', balacc_overall, ...
        'acc_std', acc_std, ...
        'balacc_std', balacc_std, ...
        'recall_class', recall_class, ...
        'cm', {cm}, ...
        'metrics', {metrics}, ...
        'chance_bal', chance_bal, ...
        'classes', classes, ...
        'increment_times_aligned', increment_times);
end