clearvars
close all
clc

session = "417"; 
v_trials = 1:32;

filename_excel = "C:\Users\manue\OneDrive - Scuola Superiore Sant'Anna\Documents\PhD - Manuela Uliano\12_Period abroad\02_Chalmers\05_Experiment\03_Data_gaze\BCI03\Notes.xlsx";
opts = detectImportOptions(filename_excel);
opts = setvartype(opts, opts.VariableNames, "string");
T = readtable(filename_excel, opts);

conditions = ["Gaze-on-target"];
    % "Free-gaze", "Gaze-on-center", "Gaze-on-target", "Gaze-only"];
filename_neural = {'../00_Data_extraction/BCI03_Session_417/controlled_BCI03.mat'};

for condition = 1:numel(conditions)
    S = load(filename_neural{condition});
    data = S.data;
    rows = (T.Condition == conditions(condition) & T.Session == session);
    for k = 1:numel(data)
        for j = 1:numel(data(k).Data)
            nO = numel(data(k).Data(j).Original);
            nI = numel(data(k).Data(j).Interp);
    
            if nO > 0
                [data(k).Data(j).Original(1:nO).Excluded] = deal(0);
            end
            if nI > 0
                [data(k).Data(j).Interp(1:nI).Excluded] = deal(0);
            end
        end
    end
    
    sets_in_table = unique(T.Set(rows));
    for s = 1:numel(sets_in_table)
        setVal = sets_in_table(s);
    
        idx = find(str2double(string([data.Set])) == str2double(setVal), 1);
        if isempty(idx)
            warning("Set %s presente in Excel ma non trovato in data.", string(setVal));
            continue
        end
    
        theseRows = find(rows & (T.Set == setVal));
        for r = theseRows'
            goodStr = strtrim(T.GoodTrials(r));
            if goodStr == "all" 
                continue
            end

            if ismissing(goodStr)
                for j = 1:numel(data(idx).Data)
                    [data(idx).Data(j).Original(:).Excluded] = deal(1);
                    [data(idx).Data(j).Interp(:).Excluded] = deal(1);
                end
                continue
            end
    
            excl = setdiff(v_trials, str2double(strtrim(split(goodStr, ","))));
            excl = excl(~isnan(excl));          % rimuovi NaN
            excl = unique(excl(:))';            % unici
    
            for j = 1:numel(data(idx).Data)
                nO = numel(data(idx).Data(j).Original);
                nI = numel(data(idx).Data(j).Interp);
    
                exclO = excl(excl >= 1 & excl <= nO);
                exclI = excl(excl >= 1 & excl <= nI);
    
                if ~isempty(exclO)
                    [data(idx).Data(j).Original(exclO).Excluded] = deal(1);
                end
                if ~isempty(exclI)
                    [data(idx).Data(j).Interp(exclI).Excluded] = deal(1);
                end
            end
        end
    end

    [folder, name, ext] = fileparts(filename_neural{condition});
    outputFile = fullfile(folder, name + "_exclUpdated" + ext);
    save(outputFile, "data", "-v7.3");
end 