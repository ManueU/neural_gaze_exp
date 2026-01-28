clearvars
close all
clc

filename_excel = "C:\Users\manue\OneDrive - Scuola Superiore Sant'Anna\Documents\PhD - Manuela Uliano\12_Period abroad\02_Chalmers\05_Experiment\03_Data_gaze\BCI02\00_Session_783\Notes.xlsx";
T = readtable(filename_excel, "TextType","string");

conditions = ["Free-gaze", "Gaze-on-center", "Gaze-on-target", "Gaze-only"];
filename_neural = {'../00_Data_extraction/free-gaze_BCI02_withtracker.mat',... 
                   '../00_Data_extraction/motor_BCI02_withtracker.mat',...
                   '../00_Data_extraction/controlled_BCI02_withtracker.mat',...
                   '../00_Data_extraction/gaze_BCI02_withtracker.mat'};

for condition = 1:numel(conditions)
    S = load(filename_neural{condition});
    data = S.data;
    rows = (T.Condition == conditions(condition));
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
    
        idx = find(str2double(string([data.Set])) == setVal, 1);
        if isempty(idx)
            warning("Set %s presente in Excel ma non trovato in data.", string(setVal));
            continue
        end
    
        theseRows = find(rows & (T.Set == setVal));
        for r = theseRows'
            excludedStr = strtrim(T.ExcludedTrials(r));
            if excludedStr == "" || ismissing(excludedStr)
                continue
            end
    
            excl = str2double(strtrim(split(excludedStr, ",")));
            excl = excl(~isnan(excl));          % rimuovi NaN
            excl = unique(excl(:))';            % unici
    
            if isempty(excl)
                continue
            end
    
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