function [cleanFactors,cleanKins,cleanMotor,cleanIDX] = CleanData(obj,TS,bad_epochs,varargin)
if isfield(TS{1}{1},'All')
    act_field = 'All';
else
    act_field = 'Motor';
end

if nargin < 3; bad_epochs = []; end
if ischar(bad_epochs)
    badeps = {bad_epochs};
elseif iscell(bad_epochs)
    badeps = bad_epochs;
else
    badeps = [];
end

kindims = obj.KinematicDims(obj.KinematicDims~=7);

epnames = TS{1}{1}.epochs.names;
badeploc = false(length(badeps),length(epnames));
for i = 1:length(badeps)
    badeploc(i,:) = cellfun(@(x) contains(x,badeps{i},'IgnoreCase',true),epnames);
end
badeploc = any(badeploc,1);

fac_cell = cellfun(@(x) x.epochs.Factors,TS{1}','uni',0); fac_cellmat = vertcat(fac_cell{:});
mot_cell = cellfun(@(x) x.epochs.(act_field),TS{1}','uni',0); mot_cellmat = vertcat(mot_cell{:});
kin_cell = cellfun(@(x) x.epochs.Kin.Vel,TS{1}','uni',0); kin_cellmat = vertcat(kin_cell{:});
kin_cellmat = cellfun(@(x) x(kindims,:),kin_cellmat,'uni',0);

preseps = find(cellfun(@(x) contains(x,'pres','IgnoreCase',true),epnames));
idx_cellmat = cellfun(@(x) true(1,size(x,2)),fac_cellmat,'uni',0);
for i = 1:length(preseps) % presentation epochs
    ei = preseps(i);
    for j = 1:size(kin_cellmat,1) % trials
        no_kin = sum(kin_cellmat{j,ei}.^2,1)<eps; % Find presentation-phase times with no kinematics
        
        if ~isempty(kin_cellmat{j,ei})
            kin_cellmat{j,ei} = kin_cellmat{j,ei}(:,no_kin);
            fac_cellmat{j,ei} = fac_cellmat{j,ei}(:,no_kin);
            mot_cellmat{j,ei} = mot_cellmat{j,ei}(:,no_kin);
            idx_cellmat{j,ei}(~no_kin) = false;
        end
    end
end
fiteps = ~badeploc;

cleanFactors = cell2mat(reshape(fac_cellmat(:,fiteps)',1,[]))';
cleanKins = cell2mat(reshape(kin_cellmat(:,fiteps)',1,[]))';
cleanMotor = cell2mat(reshape(mot_cellmat(:,fiteps)',1,[]))';

badfitidx = any(isnan([cleanFactors,cleanKins,cleanMotor]),2);

idx_cellmat(:,~fiteps) = cellfun(@(x) false(1,size(x,2)),idx_cellmat(:,~fiteps),'uni',0);
cleanIDX = cell2mat(reshape(idx_cellmat',1,[]))';

cleanFactors(badfitidx,:) = [];
cleanKins(badfitidx,:) = [];
cleanMotor(badfitidx,:) = [];
cleanIDX(badfitidx,:) = false;

end

