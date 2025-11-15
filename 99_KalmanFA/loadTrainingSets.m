function [obj,Data,ValData] = loadTrainingSets(obj,data_files,opt)
% [obj,Data, ValData] = Decoders.loadTrainingSets(obj, data_files)
% Pick a validation set and load training sets

if ~exist('data_files','var') || isempty(data_files)
    dataPath = fullfile(getClimberDataPath());
    [files, pathname]=uigetfile([dataPath,filesep,'*.bin'],'Please Select Files to load','MultiSelect','on');
    data_files = cellfun(@(x) fullfile(pathname,x),files,'UniformOutput',false);
end
if ~exist('opt','var') || isempty(opt)
    opt = struct;
end

% pick validation set
default_val_ratio = 0.2;
use_for_val = false(size(data_files));
if isfield(opt,'config') && isfield(opt.config.task_state_config,'Validation_Set')
    use_for_val = false(size(data_files));
    use_for_val(opt.config.task_state_config.Validation_Set) = true;  %%%%%NOTE: This will give unexpected results if trials are so long that they are split between multiple files!
elseif isfield(opt,'use_validation_set') && opt.use_validation_set %pull out last block of each set to use for validation
    use_for_val = false(size(data_files));
    o = regexp(data_files,'Set(?<set>[0-9]*)\.Block(?<block>[0-9]*)','once','tokens');
    sets = cellfun(@(x) str2double(x{1}),o);
    blocks = cellfun(@(x) str2double(x{2}),o);
    
    usets = unique(sets);
    for a = 1:length(usets)
        last_block_per_set = max(blocks(sets==usets(a)));
        use_for_val(sets==usets(a) & blocks == last_block_per_set)=true;
    end
elseif isfield(opt,'use_random_validation_set') && opt.use_random_validation_set %select random validation set
    t = randperm(length(use_for_val));
    use_for_val(t(1:floor(length(use_for_val)*default_val_ratio))) = true;
elseif isfield(opt, 'use_contiguous_validation_set') && opt.use_contiguous_validation_set
    tail_length = round(length(data_files) * default_val_ratio);
    use_for_val(:, end - tail_length + 1: end) = true;
end

% load training data
fprintf('Loading Saved Files and Converting...')
Data = prepData('files',data_files(~use_for_val));

% load validation data
if sum(use_for_val)>0
    ValData = prepData('files',data_files(use_for_val));
else
    ValData = [];
end

% save .mat version of training and validation data sets
if isfield(opt,'training_data_file') && ~isempty(opt.training_data_file)
    save(opt.training_data_file,'Data','ValData', '-v7.3')
end
fprintf('done\n')

end

