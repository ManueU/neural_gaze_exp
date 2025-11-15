function [obj, weights, R2] = Train(obj, data_files, opt)
global RTMA
% [OLE, weights] = Train( OLE, data_files, opt)
% Trains Decoders.RIOLEd decoder from data_files. Options can be obtained from XM
% struct

if ~exist('data_files','var')
    data_files = '';
end
if ~exist('opt','var')
    opt = struct;
end
if ~isfield(opt,'IncludeGoodTrialsOnly')
    opt.IncludeGoodTrialsOnly = true;
end
if ~isfield(opt,'use_random_validation_set')
    opt.use_random_validation_set = true;
end
if ~isfield(opt, 'report_validation_metrics') 
    % By default, CLIMBER code has been reporting train/val contaminated metrics
    % This change isn't a sophisticated check for inner-contamination (e.g. consistent cross-val or test set)
    % Just a bare minimum that we report on held-out data
    opt.report_validation_metrics = false;
end
if ~isfield(opt, 'report_r2_coeff_of_det')
    opt.report_r2_coeff_of_det = false;
end

% load training and validation data
[obj,Data, ValData] = obj.loadTrainingSets(data_files, opt);

% prep for test
if ~isempty(ValData)
    [Measurement_val, Kinematics_val, TaskStateMasks_val, idx_val]= obj.PrepDataForTest(ValData);
    %clear ValData;
    ValData = 0; % Remove uneeded data, but do not clear or set to empty to differeniate from no validation set
end

% prep for train
obj.MeasurementMask = []; 
if any(obj.IgnoreKinematicDims) %isfield(XM,'IgnoreKinematicDims') %TODO update to positive-logic (i.e. select dims we ARE training)
    [obj,Measurement, Kinematics, TaskStateMasks,idx] = obj.PrepDataForTrain(Data,'success_only',opt.IncludeGoodTrialsOnly==1,'IgnoreKinematicDims',obj.IgnoreKinematicDims);
else
    [obj,Measurement, Kinematics, TaskStateMasks,idx] = obj.PrepDataForTrain(Data,'success_only',opt.IncludeGoodTrialsOnly==1);
end

clear Data;

% moved to Decoders.Decoder.PrepDataForTrain
% obj = SetConfig(obj,'featureMeans',featureMeans);
% obj = SetConfig(obj,'featureStds',featureStds);
if isfield(opt,'use_FR_Norm') && opt.use_FR_Norm
    obj = SetConfig(obj,'summedFR',nanmean(sum(Measurement,2)));
end

% handle force calibration
if isfield(opt,'config') && isfield(opt.config.task_state_config,'force_calib') &&   ~isempty(find(not(cellfun('isempty',strfind(TaskStateMasks.states,'ForcePres'))), 1))  %for inserting intended force kinematics for decoding
    obj = SetConfig(obj,'force_decode',true);
    forceStateNum = find(not(cellfun('isempty',strfind(TaskStateMasks.states,'ForcePres')))); %Finds any state that starts with "ForcePres"
    
    intendedForceLevel = zeros(1,length(TaskStateMasks.force_calib));
    if length(forceStateNum) > 1
        %place holder for if we use multiple ForcePres phases
    else
        forceTargSet = unique(TaskStateMasks.targ_set(TaskStateMasks.state_num == forceStateNum)); %Finds the target sets used during ForcePres states
        startForcePres = [0 int32(TaskStateMasks.state_num(2:end)==forceStateNum)-int32(TaskStateMasks.state_num(1:end-1)==forceStateNum)];
        startForcePresIdxs = find(startForcePres == 1);
        for FPidx = 1:length(startForcePresIdxs)-1
            intendedForceLevel(startForcePresIdxs(FPidx):startForcePresIdxs(FPidx+1)-1) = TaskStateMasks.force_calib(startForcePresIdxs(FPidx):startForcePresIdxs(FPidx+1)-1)*...
                opt.config.task_state_config.target_configurations.(sprintf('target_set_%d',forceTargSet)).force{TaskStateMasks.targ_idx(startForcePresIdxs(FPidx))};
        end
        intendedForceLevel(startForcePresIdxs(end):end) = TaskStateMasks.force_calib(startForcePresIdxs(end):end)*...
            opt.config.task_state_config.target_configurations.(sprintf('target_set_%d',forceTargSet)).force{TaskStateMasks.targ_idx(startForcePresIdxs(end))};
        forceDim = RTMA.defines.GRIP_DIMS_R+7;
        Kinematics(:,forceDim) = intendedForceLevel(~idx)';      %Make intended force level kinematic dimension 8 (first unused dimension for 1D grasp case)
    end
    
    if ~isempty(ValData)
        intendedForceLevel = zeros(1,length(TaskStateMasks_val.force_calib));
        if length(forceStateNum) > 1
            %place holder for if we use multiple ForcePres phases
        else
            forceTargSet = unique(TaskStateMasks_val.targ_set(TaskStateMasks_val.state_num == forceStateNum)); %Finds the target sets used during ForcePres states
            startForcePres = [0 int32(TaskStateMasks_val.state_num(2:end)==forceStateNum)-int32(TaskStateMasks_val.state_num(1:end-1)==forceStateNum)];
            startForcePresIdxs = find(startForcePres == 1);
            for FPidx = 1:length(startForcePresIdxs)-1
                intendedForceLevel(startForcePresIdxs(FPidx):startForcePresIdxs(FPidx+1)-1) = TaskStateMasks_val.force_calib(startForcePresIdxs(FPidx):startForcePresIdxs(FPidx+1)-1)*...
                    opt.config.task_state_config.target_configurations.(sprintf('target_set_%d',forceTargSet)).force{TaskStateMasks_val.targ_idx(startForcePresIdxs(FPidx))};
            end
            intendedForceLevel(startForcePresIdxs(end):end) = TaskStateMasks_val.force_calib(startForcePresIdxs(end):end)*...
                opt.config.task_state_config.target_configurations.(sprintf('target_set_%d',forceTargSet)).force{TaskStateMasks_val.targ_idx(startForcePresIdxs(end))};
            Kinematics_val(:,8) = intendedForceLevel(~idx_val)';      %Make intended force level kinematic dimension 8 (first unused dimension for 1D grasp case)
        end
        
    end
elseif isfield(opt,'config') && ((isfield(opt.config.task_state_config,'force_calib') && opt.config.task_state_config.force_calib == 1) || (isfield(opt.config.task_state_config,'gripperTask') && opt.config.task_state_config.gripperTask==1))
    obj = SetConfig(obj,'force_decode',true);
end

fprintf('Using %d features\n',sum(obj.MeasurementMask))
%%% TRAINING:
fprintf('Calling Decoder Train Method... ')
if ~isempty(ValData)
    [obj, weights] = Train_helper(obj,Measurement,Kinematics,Measurement_val,Kinematics_val, opt);
    % Report accuracy:
    if opt.report_validation_metrics
        [~,R2] = obj.SimpleValidate([Measurement_val],[Kinematics_val], false, opt.report_r2_coeff_of_det);
    else
        [~,R2] = obj.SimpleValidate([Measurement; Measurement_val],[Kinematics; Kinematics_val], false, opt.report_r2_coeff_of_det);
    end
else
    [obj, weights] = Train_helper(obj,Measurement,Kinematics, [], [], opt);
    [~,R2] = obj.SimpleValidate(Measurement,Kinematics, false, opt.report_r2_coeff_of_det);
end

fprintf('Training Complete.\n')
end


function [obj, weights] = Train_helper( obj, NeuralSig, KinematicSig, NeuralSig_val, KinematicSig_val, opt)
global RTMA
% [OLE, weights] = Train( OLE, NeuralSig,KinematicSig)
%

% Brian Wodlinger based on Meel Velliste 9/19/2009

if isempty(RTMA)
    RTMA = LoadRtmaConfig;
end

obj.Kmeans = nanmean(KinematicSig,1);
obj.Kstd = nanstd(KinematicSig,[],1);

IgnoreIdx = obj.IgnoreKinematicDims;
IgnoreIdx(size(KinematicSig,2)+1:end) = []; %limit to same size as Kinematic sig, since we are probably not using the last 8 dims.
IgnoreIdx = IgnoreIdx | obj.Kstd==0 | isnan(obj.Kstd);
obj.Kmeans(IgnoreIdx==1) = 0;
obj.Kstd(IgnoreIdx==1) = 1;

KinematicSig = KinematicSig -  ones(size(KinematicSig,1),1)*obj.Kmeans;
KinematicSig = KinematicSig./(ones(size(KinematicSig,1),1)*obj.Kstd);
KinematicSig_val = KinematicSig_val -  ones(size(KinematicSig_val,1),1)*obj.Kmeans;
KinematicSig_val = KinematicSig_val./(ones(size(KinematicSig_val,1),1)*obj.Kstd);



%standardize
if isempty(obj.featureMeans)
    obj.featureMeans = nanmean(NeuralSig,1);
    obj.featureStds = nanstd(NeuralSig,[],1);
end
NeuralSig = NeuralSig -  ones(size(NeuralSig,1),1)*obj.featureMeans;
NeuralSig = NeuralSig./(ones(size(NeuralSig,1),1)*obj.featureStds);

if isempty(obj.MeasurementMask)
    obj.MeasurementMask = ones(1,size(NeuralSig,2));
end
NeuralSig = NeuralSig(:,obj.MeasurementMask==1);  % Features to reject before we even do any training

%standardize validation set:
if exist('NeuralSig_val','var') && ~isempty(NeuralSig_val)
    do_validation = true;
    NeuralSig_val = NeuralSig_val -  ones(size(NeuralSig_val,1),1)*obj.featureMeans;
    NeuralSig_val = NeuralSig_val./(ones(size(NeuralSig_val,1),1)*obj.featureStds);
    
    NeuralSig_val = NeuralSig_val(:,obj.MeasurementMask==1);  % Features to reject before we even do any training
else
    do_validation = false;
end

%choose training and validation sets:
num_obs = size(NeuralSig,1); %after duplicating for balance
K_train = KinematicSig;
F_train = NeuralSig;
if do_validation
    K_val = KinematicSig_val;
    F_val = NeuralSig_val;
else
    K_val = KinematicSig;
    F_val = NeuralSig;
end
invSigma = zeros(size(F_train,2),size(F_train,2),RTMA.defines.NUM_DOMAINS);
Kt = ControlSpaceToDomainCells(K_train);
Kv = ControlSpaceToDomainCells(K_val);
K = ControlSpaceToDomainCells([KinematicSig; KinematicSig_val]);
I = ControlSpaceToDomainCells(IgnoreIdx);
K_train = K_train(:,~IgnoreIdx);
K_val = K_val(:,~IgnoreIdx);


%calculate sigma matrix:
for i = 1:size(F_train,2)
    for d = 1:RTMA.defines.NUM_DOMAINS
        [~, ~, Residuals, ~, ~] = regress(single(F_train(:,i)), [ones(size(Kt{d},1), 1) Kt{d}(:,~I{d})]);
        Sigma = var(Residuals);
        invSigma(i,i,d) = 1/Sigma;
    end
end

%determine regularization coefficient(s) to use:
lambda = [0 logspace(-3,6,100)];
lambda1= [0 logspace(-3,6,100)];
idx1=1:length(lambda1);
idx = 1:length(lambda);
metric = inf(max(idx),max(idx1));
fprintf('\n')
counter('%d/%d',1,max(idx))
s = warning('off');
if ~do_validation
    idx = 1;
    idx1=1;
    lambda = 0;
    lambda1=0;
end

for b = idx1 %limit to a features
    if b>1 && (mod(b,10)==1 || b==max(idx1))
        counter('%d/%d',b,max(idx1))
    end
    W1 = ([ones(size(K_train,1), 1) K_train]'*[ones(size(K_train,1), 1) K_train] + lambda1(b)*eye(size(K_train,2)+1))\[ones(size(K_train,1), 1) K_train]'*F_train;
    baselines = W1(1,:);
    W1 = W1(2:end,:);
    repbaselines = repmat(baselines,[size(F_val,1) 1]); %for speed
    F_val_minus_baselines = (F_val-repbaselines);
    
    W1_ = zeros(size(W1,2),length(IgnoreIdx));
    W1_(:,~IgnoreIdx) = W1';
    W1 = ControlSpaceToDomainCells(W1_);
    W1 = cellfun(@(x,i) x(:,~i),W1,I,'UniformOutput',false);
    domains_to_use = find(~cellfun(@isempty,W1));
    w1invSigma = cell(1,max(domains_to_use));
    for d = domains_to_use
        w1invSigma{d} = W1{d}'*invSigma(:,:,d);
    end
    
    for a = idx
        W = [];
        for d = domains_to_use
            Wt = (pinv(w1invSigma{d}*W1{d} + lambda(a)*eye(size(W1{d}',1)))*w1invSigma{d})';
            W = [W Wt];
        end
        
        m = corrcoef(F_val_minus_baselines*W,K_val);
        metric(a,b) = m(1,2).^2; %norm(F_val*W-K_val);
    end
    
end

if sum(~isinf(metric))>1
    figure;semilogx(lambda,metric,'*-');
    xlabel('lambda'); ylabel('R2'); 
end

warning(s);

[v, reg_parm_to_use]=max(metric(:)); %lowest error / highest r2
[reg_parm_to_use, reg_parm_to_use1] = ind2sub(size(metric),reg_parm_to_use);

fprintf('Using lambda = %f, lambda1 = %f for %d observations... \n',lambda(reg_parm_to_use),lambda1(reg_parm_to_use1),size(NeuralSig,1))
lambda = lambda(reg_parm_to_use);
lambda1 = lambda1(reg_parm_to_use1);

%final train with full dataset:
if do_validation
    if ~opt.report_validation_metrics
        NeuralSig = [NeuralSig; NeuralSig_val];
        KinematicSig = [KinematicSig; KinematicSig_val];
    end
end
KinematicSig = KinematicSig(:,~IgnoreIdx);

W = ([ones(size(KinematicSig,1), 1) KinematicSig]'*[ones(size(KinematicSig,1), 1) KinematicSig] + lambda1*eye(size(KinematicSig,2)+1))\[ones(size(KinematicSig,1), 1) KinematicSig]'*NeuralSig;
baselines = W(1,:);
W1 = W(2:end,:);
W = [];
W1_ = zeros(size(W1,2),length(IgnoreIdx));
W1_(:,~IgnoreIdx) = W1';
W1 = ControlSpaceToDomainCells(W1_);
W1 = cellfun(@(x,i) x(:,~i),W1,I,'UniformOutput',false);
domains_to_use = find(~cellfun(@isempty,W1));
for d = domains_to_use
    w1invSigma = W1{d}'*invSigma(:,:,d);
    Wt = (pinv(w1invSigma*W1{d} + lambda*eye(size(W1{d}',1)))*w1invSigma)';
    W = [W Wt];
end

weights = W;

%rescale weights so gain is reasonable:
scale = sqrt(mean((KinematicSig).^2,1)) ./ sqrt(mean(((NeuralSig-repmat(baselines,[size(NeuralSig,1) 1]))*weights).^2,1));
weights = weights .* (ones(size(weights,1),1)*scale);

obj.baselines = baselines;
obj.weights = weights;

%Report accuracy:
if opt.report_validation_metrics
    KinematicSig_val = KinematicSig_val(:,~IgnoreIdx);
    out = (NeuralSig_val-repmat(baselines,[size(NeuralSig_val,1) 1]))*weights;
    r = corrcoef(out,KinematicSig_val);
    normerr = norm(out-KinematicSig_val);
else
    out = (NeuralSig-repmat(baselines,[size(NeuralSig,1) 1]))*weights;
    r = corrcoef(out,KinematicSig);
    normerr = norm(out-KinematicSig);
end

%scale back up if we have ignored any kinematic dims:
obj.weights = zeros(size(weights,1),length(IgnoreIdx));
obj.weights(:,~(IgnoreIdx==1)) = weights;
end
