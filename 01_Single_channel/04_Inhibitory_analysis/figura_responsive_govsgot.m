% Richiede nel workspace:
%   T_go   : tabella condizione GO
%   T_got  : tabella condizione GOT
%
% Variabili richieste nelle tabelle:
%   target_id
%   array_name
%   channel_global   (oppure channel_local)

clearvars -except T_go T_got

%% --- Controlli base ---
requiredVars = {'target_id','array_name'};
for i = 1:numel(requiredVars)
    assert(ismember(requiredVars{i}, T_go.Properties.VariableNames), ...
        'T_go non contiene la variabile "%s".', requiredVars{i});
    assert(ismember(requiredVars{i}, T_got.Properties.VariableNames), ...
        'T_got non contiene la variabile "%s".', requiredVars{i});
end

if ismember('channel_global', T_go.Properties.VariableNames) && ...
   ismember('channel_global', T_got.Properties.VariableNames)
    channelVar = 'channel_global';
elseif ismember('channel_local', T_go.Properties.VariableNames) && ...
       ismember('channel_local', T_got.Properties.VariableNames)
    channelVar = 'channel_local';
else
    error('Nessuna variabile canale trovata: né channel_global né channel_local.');
end

%% --- Uniforma tipi ---
T_go.array_name  = string(T_go.array_name);
T_got.array_name = string(T_got.array_name);

if iscell(T_go.target_id),  T_go.target_id  = string(T_go.target_id);  end
if iscell(T_got.target_id), T_got.target_id = string(T_got.target_id); end

%% --- Lista target ---
targets_go  = unique(T_go.target_id);
targets_got = unique(T_got.target_id);
targets     = union(targets_go, targets_got, 'stable');

nTargets = numel(targets);

n_common   = zeros(nTargets,1);
n_only_go  = zeros(nTargets,1);
n_only_got = zeros(nTargets,1);

%% --- Conta overlap per target ---
for i = 1:nTargets
    tgt = targets(i);

    go_t  = T_go(T_go.target_id == tgt, :);
    got_t = T_got(T_got.target_id == tgt, :);

    key_go  = unique(go_t.array_name  + "_" + string(go_t.(channelVar)));
    key_got = unique(got_t.array_name + "_" + string(got_t.(channelVar)));

    common_keys   = intersect(key_go, key_got);
    only_go_keys  = setdiff(key_go, key_got);
    only_got_keys = setdiff(key_got, key_go);

    n_common(i)   = numel(common_keys);
    n_only_go(i)  = numel(only_go_keys);
    n_only_got(i) = numel(only_got_keys);
end

%% --- Matrice stacked bar ---
Y = [n_common, n_only_go, n_only_got];

%% --- Plot ---
figure('Color','w','Name','Responsive channel overlap by target');
b = bar(Y, 'stacked', 'LineWidth', 1.0);

% Colori semplici e leggibili
b(1).FaceColor = [0.55 0.55 0.55];   % common = grigio
b(2).FaceColor = [0.20 0.45 0.85];   % only GO = blu
b(3).FaceColor = [0.85 0.40 0.20];   % only GOT = arancio

ax = gca;
ax.XTick = 1:nTargets;
ax.XTickLabel = cellstr("T" + targets);
ax.Box = 'off';
ax.FontSize = 11;

xlabel('Target')
ylabel('Responsive channels')
title('Responsive channels across conditions')

legend({'Common','Only GO','Only GOT'}, ...
    'Location','northoutside', ...
    'Orientation','horizontal', ...
    'Box','off')

grid on
ax.YGrid = 'on';
ax.XGrid = 'off';

%% --- Etichetta del totale sopra ogni barra ---
totals = sum(Y, 2);
hold on
for i = 1:nTargets
    text(i, totals(i) + 0.3, sprintf('n=%d', totals(i)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 10);
end
hold off