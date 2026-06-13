%% Selezione dei canali comuni tra due condizioni, separatamente per target
% File:
%   1) inhibitory_analysis_go.mat
%   2) inhibitory_analysis_got.mat
%
% Output:
% - common_by_target: struct con, per ogni target,
%       .go   = righe della condizione GO con canali comuni
%       .got  = righe della condizione GOT con canali comuni

clear
close all
clc

%% --- Load tables from .mat files ---
goMat  = load('inhibitory_analysis_go.mat');
gotMat = load('inhibitory_analysis_got.mat');

T_go  = extractTableFromMat(goMat);
T_got = extractTableFromMat(gotMat);

channelVar = 'channel_global'; 

%% --- Targets presenti almeno in una delle due condizioni ---
targets_go  = unique(T_go.target_id);
targets_got = unique(T_got.target_id);
all_targets = union(targets_go, targets_got);

%% --- Select common channels for each target ---
common_by_target = struct();

for iT = 1:numel(all_targets)
    tgt = all_targets(iT);

    % Righe del target corrente
    go_t  = T_go(T_go.target_id == tgt, :);
    got_t = T_got(T_got.target_id == tgt, :);

    % Chiave univoca canale = array_name + channel
    key_go  = string(go_t.(channelVar));
    key_got = string(got_t.(channelVar));

    % Intersezione dei canali tra le due condizioni
    common_keys = intersect(key_go, key_got);

    % Mantieni solo i canali comuni
    go_common  = go_t(ismember(key_go, common_keys), :);
    got_common = got_t(ismember(key_got, common_keys), :);

    % Ordina per canale per avere tabelle allineate
    go_common  = sortrows(go_common);
    got_common = sortrows(got_common);

    % Nome campo valido per struct
    fieldName = matlab.lang.makeValidName("target_" + string(tgt));

    common_by_target.(fieldName).target_id = tgt;
    common_by_target.(fieldName).go  = go_common;
    common_by_target.(fieldName).got = got_common;
    common_by_target.(fieldName).n_common_channels = numel(common_keys);
end

%% --- Display summary ---
disp('Numero di canali comuni per target:');
fn = fieldnames(common_by_target);
for i = 1:numel(fn)
    fprintf('%s -> %d canali comuni\n', ...
        string(common_by_target.(fn{i}).target_id), ...
        common_by_target.(fn{i}).n_common_channels);
end


%% Boxplot per ciascun target usando solo i canali comuni
% Richiede nel workspace:
%   common_by_target
%
% Ogni pannello mostra, per un target:
%   - boxplot condizione GO
%   - boxplot condizione GOT
% costruiti solo sui canali presenti in entrambe le condizioni.

measureVar = 'peak_inh_amp';

fn = fieldnames(common_by_target);
nTargets = numel(fn);

figure('Color','w');
hold on

data = {};
pos = [];
targetCenters = zeros(nTargets,1);
targetLabels = cell(nTargets,1);

for i = 1:nTargets
    S = common_by_target.(fn{i});

    goVals  = S.go.(measureVar);
    gotVals = S.got.(measureVar);

    goVals  = goVals(~isnan(goVals));
    gotVals = gotVals(~isnan(gotVals));

    % due posizioni per ogni target
    p1 = (i-1)*3 + 1;   % GO
    p2 = (i-1)*3 + 2;   % GOT

    if ~isempty(goVals)
        data{end+1} = goVals;
        pos(end+1) = p1;
    end

    if ~isempty(gotVals)
        data{end+1} = gotVals;
        pos(end+1) = p2;
    end

    targetCenters(i) = mean([p1 p2]);
    targetLabels{i} = sprintf('Target %s', string(S.target_id));
end

% boxplot richiede vettore unico
allY = [];
allPos = [];
for k = 1:numel(data)
    allY = [allY; data{k}(:)];
    allPos = [allPos; repmat(pos(k), numel(data{k}), 1)];
end

boxplot(allY, allPos, 'Positions', unique(allPos), 'Widths', 0.6)

xlim([0 max(pos)+1])
set(gca, 'XTick', targetCenters, 'XTickLabel', targetLabels)

ylabel('Peak amplitude relative to baseline (Hz)')
title('Gaze-only vs Gaze-on-target')
grid on

% etichette GO/GOT sotto ciascun box
yl = ylim;
for i = 1:nTargets
    p1 = (i-1)*3 + 1;
    p2 = (i-1)*3 + 2;

    text(p1, yl(1) - 0.05*(yl(2)-yl(1)), 'GO', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
    text(p2, yl(1) - 0.05*(yl(2)-yl(1)), 'GOT', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
end
box off
hold off


%% --- Local function ---
function T = extractTableFromMat(S)
    % Estrae la prima variabile di classe table trovata nel file .mat
    vars = fieldnames(S);
    T = [];
    for k = 1:numel(vars)
        if istable(S.(vars{k}))
            T = S.(vars{k});
            return;
        end
    end
    error('Nel file .mat non è stata trovata alcuna tabella.');
end