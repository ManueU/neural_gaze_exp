close all
clear samplesAfter
clear samplesBefore

states = set04.TaskStateMasks.state_name;   % cell array
numericStates = zeros(size(states));        % preallocazione

for i = 1:numel(states)
    switch states{i}   
        case 'Center'
            numericStates(i) = 1;
        case 'Pres12'
            numericStates(i) = 2;
        case 'Reach'
            numericStates(i) = 3;      
        case 'Remove'
            numericStates(i) = 10;     
        otherwise
            numericStates(i) = NaN;
    end
end

trialNum = set04.trial_num;  
d = diff(trialNum);

% Gli indici dove cambia sono quelli in cui diff ≠ 0
changeIdx = find(d ~= 0) + 1; 


figure()
plot(numericStates, 'b'), hold on
xline(changeIdx, 'r')


%%

state = numericStates;      % vettore degli stati (1..4) per ogni campione
redIdx = changeIdx;       % vettore con l’indice (in campioni) di ogni linea rossa

N = numel(state);

% Indici dove la traccia blu cambia livello (gradini)
stepIdx = find(diff(state) ~= 0) + 1;   % inizio di ogni nuovo livello

nR = numel(redIdx);
samplesBefore = nan(nR,1);
samplesAfter  = nan(nR,1);

for k = 1:nR
    i = redIdx(k);   % campione corrispondente alla k-esima linea rossa

    % ---- gradino precedente ----
    idxPrev = find(stepIdx < i, 1, 'last');
    if ~isempty(idxPrev)
        prevStep = stepIdx(idxPrev);
        samplesBefore(k) = i - prevStep;      % campioni tra gradino precedente e linea rossa
    end

    % ---- gradino successivo ----
    idxNext = find(stepIdx > i, 1, 'first');
    if ~isempty(idxNext)
        nextStep = stepIdx(idxNext);
        samplesAfter(k) = nextStep - i;       % campioni tra linea rossa e gradino successivo
    end
end
