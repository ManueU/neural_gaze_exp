%% Population-level boxplots: GO vs GOT per ciascun target
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

T_go.array_name  = string(T_go.array_name);
T_got.array_name = string(T_got.array_name);

if iscell(T_go.target_id),  T_go.target_id  = string(T_go.target_id);  end
if iscell(T_got.target_id), T_got.target_id = string(T_got.target_id); end

targets = union(unique(T_go.target_id), unique(T_got.target_id), 'stable');
nTargets = numel(targets);

%% ---------------------------
% 2) Costruzione dati per boxplot
% ---------------------------
allData   = [];
groupNum  = [];
positions = [];
labelsBox = {};
meansBox  = [];
nsBox     = [];

% posizioni: per ogni target due box
% es. T1 -> 1,2 ; T2 -> 4,5 ; T3 -> 7,8 ...
gap = 1;
posGO  = zeros(nTargets,1);
posGOT = zeros(nTargets,1);
targetCenters = zeros(nTargets,1);
targetLabels = cell(nTargets,1);

currentPos = 1;

for i = 1:nTargets
    tgt = targets(i);

    goVals  = T_go{T_go.target_id == tgt, measureVar};
    gotVals = T_got{T_got.target_id == tgt, measureVar};

    goVals  = goVals(~isnan(goVals));
    gotVals = gotVals(~isnan(gotVals));

    posGO(i)  = currentPos;
    posGOT(i) = currentPos + 1;
    targetCenters(i) = mean([posGO(i), posGOT(i)]);
    targetLabels{i} = sprintf('%s', string(tgt));

    if ~isempty(goVals)
        allData   = [allData; goVals(:)];
        groupNum  = [groupNum; repmat(2*i-1, numel(goVals), 1)];
        positions = [positions; repmat(posGO(i), numel(goVals), 1)];
        labelsBox{end+1} = sprintf('GO_T%s', string(tgt));
        meansBox(end+1)  = mean(goVals);
        nsBox(end+1)     = numel(goVals);
    else
        labelsBox{end+1} = sprintf('GO_T%s', string(tgt));
        meansBox(end+1)  = NaN;
        nsBox(end+1)     = 0;
    end

    if ~isempty(gotVals)
        allData   = [allData; gotVals(:)];
        groupNum  = [groupNum; repmat(2*i, numel(gotVals), 1)];
        positions = [positions; repmat(posGOT(i), numel(gotVals), 1)];
        labelsBox{end+1} = sprintf('GOT_T%s', string(tgt));
        meansBox(end+1)  = mean(gotVals);
        nsBox(end+1)     = numel(gotVals);
    else
        labelsBox{end+1} = sprintf('GOT_T%s', string(tgt));
        meansBox(end+1)  = NaN;
        nsBox(end+1)     = 0;
    end

    currentPos = currentPos + 2 + gap;
end

boxPositions = reshape([posGO posGOT].', [], 1);

%% ---------------------------
% 3) Plot
% ---------------------------
figure('Color','w','Position',[100 100 1100 450]);
hold on

boxplot(allData, positions, ...
    'Positions', unique(positions), ...
    'Widths', 0.65, ...
    'Colors', [0 0 0], ...
    'Symbol', '', ...
    'Whisker', 1.5);

% Colori
colGO  = [0.20 0.45 0.85];
colGOT = [0.90 0.45 0.15];

% Trova box e colorale alternando GO/GOT
h = findobj(gca, 'Tag', 'Box');
h = flipud(h);

for k = 1:numel(h)
    if mod(k,2)==1
        thisCol = colGO;
    else
        thisCol = colGOT;
    end

    patch(get(h(k),'XData'), get(h(k),'YData'), thisCol, ...
        'FaceAlpha', 0.60, ...
        'EdgeColor', 'k', ...
        'LineWidth', 1.0);
end

% Ridisegna mediane e whiskers sopra le patch
set(findobj(gca,'Tag','Median'), 'Color', 'k', 'LineWidth', 1.4);
set(findobj(gca,'Tag','Whisker'),'Color','k','LineWidth',1.0);
set(findobj(gca,'Tag','Upper Whisker'),'Color','k','LineWidth',1.0);
set(findobj(gca,'Tag','Lower Whisker'),'Color','k','LineWidth',1.0);
set(findobj(gca,'Tag','Upper Adjacent Value'),'Color','k','LineWidth',1.0);
set(findobj(gca,'Tag','Lower Adjacent Value'),'Color','k','LineWidth',1.0);

%% ---------------------------
% 4) Overlay punti individuali
% ---------------------------
for i = 1:nTargets
    tgt = targets(i);

    goVals  = T_go{T_go.target_id == tgt, measureVar};
    gotVals = T_got{T_got.target_id == tgt, measureVar};

    goVals  = goVals(~isnan(goVals));
    gotVals = gotVals(~isnan(gotVals));

    if ~isempty(goVals)
        xj = posGO(i) + 0.16*(rand(size(goVals)) - 0.5);
        scatter(xj, goVals, 20, ...
            'MarkerFaceColor', colGO, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.35, ...
            'MarkerFaceAlpha', 0.45, ...
            'MarkerEdgeAlpha', 0.25);
    end

    if ~isempty(gotVals)
        xj = posGOT(i) + 0.16*(rand(size(gotVals)) - 0.5);
        scatter(xj, gotVals, 20, ...
            'MarkerFaceColor', colGOT, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.35, ...
            'MarkerFaceAlpha', 0.45, ...
            'MarkerEdgeAlpha', 0.25);
    end
end

%% ---------------------------
% 5) Media come diamanti neri
% ---------------------------
for i = 1:nTargets
    goVals  = T_go{T_go.target_id == targets(i), measureVar};
    gotVals = T_got{T_got.target_id == targets(i), measureVar};

    goVals  = goVals(~isnan(goVals));
    gotVals = gotVals(~isnan(gotVals));

    if ~isempty(goVals)
        plot(posGO(i), mean(goVals), 'kd', ...
            'MarkerFaceColor', 'k', 'MarkerSize', 6)
    end

    if ~isempty(gotVals)
        plot(posGOT(i), mean(gotVals), 'kd', ...
            'MarkerFaceColor', 'k', 'MarkerSize', 6)
    end
end

%% ---------------------------
% 6) Asse e annotazioni
% ---------------------------
set(gca, ...
    'XTick', targetCenters, ...
    'XTickLabel', targetLabels, ...
    'FontSize', 11, ...
    'LineWidth', 1.0, ...
    'Box', 'off', ...
    'TickDir', 'out');

xlabel('Target')
ylabel('Peak amplitude relative to baseline (Hz)')
title('Population-level inhibitory amplitude by target and condition')

grid on
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

xlim([min(posGO)-1, max(posGOT)+1])

% etichette GO / GOT sotto ogni coppia
% yl = ylim;
% yr = yl(2)-yl(1);
% for i = 1:nTargets
%     text(posGO(i),  yl(1)-0.06*yr, 'GO', ...
%         'HorizontalAlignment','center', ...
%         'VerticalAlignment','top', ...
%         'FontSize',10, ...
%         'Color', colGO);
% 
%     text(posGOT(i), yl(1)-0.06*yr, 'GOT', ...
%         'HorizontalAlignment','center', ...
%         'VerticalAlignment','top', ...
%         'FontSize',10, ...
%         'Color', colGOT);
% end

% n sopra ogni box
yl = ylim;
yr = yl(2)-yl(1);

for i = 1:nTargets
    goVals  = T_go{T_go.target_id == targets(i), measureVar};
    gotVals = T_got{T_got.target_id == targets(i), measureVar};

    nGO  = sum(~isnan(goVals));
    nGOT = sum(~isnan(gotVals));

    text(posGO(i),  yl(2)+0.03*yr, sprintf('n=%d', nGO), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',9, ...
        'Color', colGO);

    text(posGOT(i), yl(2)+0.03*yr, sprintf('n=%d', nGOT), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',9, ...
        'Color', colGOT);
end

ylim([yl(1)-0.10*yr, yl(2)+0.10*yr])

%% ---------------------------
% 7) Legend minimale
% ---------------------------
p1 = plot(nan, nan, 's', 'MarkerSize', 8, ...
    'MarkerFaceColor', colGO, 'MarkerEdgeColor', 'k');
p2 = plot(nan, nan, 's', 'MarkerSize', 8, ...
    'MarkerFaceColor', colGOT, 'MarkerEdgeColor', 'k');
p3 = plot(nan, nan, 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

legend([p1 p2 p3], {'GO','GOT','Mean'}, ...
    'Location', 'northoutside', ...
    'Orientation', 'horizontal', ...
    'Box', 'off');

hold off