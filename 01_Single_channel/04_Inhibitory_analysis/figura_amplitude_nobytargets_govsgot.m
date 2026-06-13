% Distribuzioni aggregate:
%   common_GO, common_GOT, only_GO, only_GOT
%
% Richiede nel workspace:
%   T_go
%   T_got

clearvars -except T_go T_got

measureVar = 'peak_inh_amp';

%% ---------------------------
% 1) Controlli e setup
% ---------------------------
requiredVars = {'target_id','array_name',measureVar};
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

T_go.array_name  = string(T_go.array_name);
T_got.array_name = string(T_got.array_name);

if iscell(T_go.target_id),  T_go.target_id  = string(T_go.target_id);  end
if iscell(T_got.target_id), T_got.target_id = string(T_got.target_id); end

targets = union(unique(T_go.target_id), unique(T_got.target_id), 'stable');

%% ---------------------------
% 2) Costruzione gruppi aggregati
% ---------------------------
common_GO  = [];
common_GOT = [];
only_GO    = [];
only_GOT   = [];

for i = 1:numel(targets)
    tgt = targets(i);

    go_t  = T_go(T_go.target_id == tgt, :);
    got_t = T_got(T_got.target_id == tgt, :);

    key_go  = go_t.array_name  + "_" + string(go_t.(channelVar));
    key_got = got_t.array_name + "_" + string(got_t.(channelVar));

    common_keys   = intersect(key_go, key_got);
    only_go_keys  = setdiff(key_go, key_got);
    only_got_keys = setdiff(key_got, key_go);

    go_common  = go_t(ismember(key_go, common_keys), :);
    got_common = got_t(ismember(key_got, common_keys), :);

    go_only  = go_t(ismember(key_go, only_go_keys), :);
    got_only = got_t(ismember(key_got, only_got_keys), :);

    common_GO  = [common_GO;  go_common.(measureVar)];
    common_GOT = [common_GOT; got_common.(measureVar)];
    only_GO    = [only_GO;    go_only.(measureVar)];
    only_GOT   = [only_GOT;   got_only.(measureVar)];
end

common_GO  = common_GO(~isnan(common_GO));
common_GOT = common_GOT(~isnan(common_GOT));
only_GO    = only_GO(~isnan(only_GO));
only_GOT   = only_GOT(~isnan(only_GOT));

dataCell = {common_GO, common_GOT, only_GO, only_GOT};
labels   = {'Common GO','Common GOT','Only GO','Only GOT'};
xpos     = 1:4;

ns = cellfun(@numel, dataCell);

%% ---------------------------
% 3) Colori
% ---------------------------
% GO = blu, GOT = arancio
% common = più saturo, only = più chiaro
col_common_go  = [0.20 0.45 0.85];
col_common_got = [0.90 0.45 0.15];
col_only_go    = [0.65 0.78 0.95];
col_only_got   = [0.97 0.75 0.60];

cols = [
    col_common_go
    col_common_got
    col_only_go
    col_only_got
];

%% ---------------------------
% 4) Vettori per boxplot
% ---------------------------
allData = [];
groupNum = [];

for k = 1:numel(dataCell)
    allData  = [allData; dataCell{k}(:)];
    groupNum = [groupNum; repmat(k, numel(dataCell{k}), 1)];
end

%% ---------------------------
% 5) Plot
% ---------------------------
figure('Color','w','Position',[100 100 700 450]);
hold on

% Boxplot pulito, senza simboli outlier automatici
boxplot(allData, groupNum, ...
    'Positions', xpos, ...
    'Widths', 0.55, ...
    'Colors', [0 0 0], ...
    'Symbol', '', ...
    'Whisker', 1.5);

% Riempi manualmente le box con patch colorate
h = findobj(gca, 'Tag', 'Box');
h = flipud(h); % ordine coerente con xpos

for k = 1:min(numel(h),4)
    patch(get(h(k),'XData'), get(h(k),'YData'), cols(k,:), ...
        'FaceAlpha', 0.65, 'EdgeColor', 'k', 'LineWidth', 1.0);
end

% Ridisegna mediane e whisker sopra le patch
set(findobj(gca,'Tag','Median'), 'Color', 'k', 'LineWidth', 1.4);
set(findobj(gca,'Tag','Whisker'),'Color','k','LineWidth',1.0);
set(findobj(gca,'Tag','Upper Whisker'),'Color','k','LineWidth',1.0);
set(findobj(gca,'Tag','Lower Whisker'),'Color','k','LineWidth',1.0);
set(findobj(gca,'Tag','Upper Adjacent Value'),'Color','k','LineWidth',1.0);
set(findobj(gca,'Tag','Lower Adjacent Value'),'Color','k','LineWidth',1.0);

% Punti individuali jitterati
for k = 1:numel(dataCell)
    xj = xpos(k) + 0.16*(rand(size(dataCell{k})) - 0.5);
    scatter(xj, dataCell{k}, 22, ...
        'MarkerFaceColor', cols(k,:), ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.35, ...
        'MarkerFaceAlpha', 0.55, ...
        'MarkerEdgeAlpha', 0.35);
end

% Media come diamante nero
means = cellfun(@mean, dataCell);
plot(xpos, means, 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 6)

%% ---------------------------
% 6) Asse e annotazioni
% ---------------------------
set(gca, ...
    'XTick', xpos, ...
    'XTickLabel', labels, ...
    'FontSize', 11, ...
    'LineWidth', 1.0, ...
    'Box', 'off', ...
    'TickDir', 'out');

xtickangle(0)

ylabel('Peak amplitude relative to baseline (Hz)')
title('Distribution of inhibitory amplitude')

grid on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% Aggiungi n sopra ogni gruppo
yl = ylim;
yr = yl(2) - yl(1);
yText = yl(2) + 0.04*yr;

for k = 1:4
    text(xpos(k), yText, sprintf('n = %d', ns(k)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 10);
end

ylim([yl(1), yl(2) + 0.10*yr])

hold off