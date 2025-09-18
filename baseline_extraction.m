close all 
clc

[file, path] = uigetfile('*.mat', 'Seleziona un file MAT');
if isequal(file,0)
    error('Nessun file selezionato.');
else
    disp('Loading data...')
    load(fullfile(path, file));
end

counts = cellfun(@numel, baselineStruct.spikes.motor);
bas_m_medial = counts(1:96)./60;
bas_m_lateral = counts(97:192)./60;