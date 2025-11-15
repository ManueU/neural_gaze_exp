function [obj,Measurement,Kinematics,TaskStateMasks,idx,reserved] = PrepDataForTrain(obj,Data,varargin)
% To be called after prepData  in order to apply the filters, remove nans et cetra to get signals ready for training

parse_varargin(varargin,'success_only',false,'remove_bad_trials',true,'IgnoreKinematicDims',[]);

reserved = []; % variable data reserved for subclasses

% define scale factors for the different dimensions.
trans_factor = 1;
ori_factor = 1;
grasp_factor = .5;

% Prep data for training (copied from SimpleAdapt)
Measurement = double(Data.SpikeCount);
if isfield(Data.Kinematics,'ActualForce')	% Needs to be adjusted if MuJoCo goes >1 DOF. This currently doesn't include translation either, and assumes we want to use force
    try
        Kinematics = [Data.Kinematics.Control(:,1:6), Data.Kinematics.CommVel, Data.Kinematics.ActualForce(:,1) zeros(length(Data.Kinematics.ActualVel),22)];
    catch
        Kinematics = [zeros(length(Data.Kinematics.ActualVel),6) Data.Kinematics.ActualVel Data.Kinematics.ActualForce(:,1) zeros(length(Data.Kinematics.ActualVel),22)];
    end
else
    Kinematics = Data.Kinematics.Control;
end
TaskStateMasks = Data.TaskStateMasks;

Kinematics(:,TransIdxs) = Kinematics(:,TransIdxs)*trans_factor;
Kinematics(:,OriIdxs) = Kinematics(:,OriIdxs)*ori_factor;
Kinematics(:,GraspIdxs) = Kinematics(:,GraspIdxs)*grasp_factor;

rtma = LoadRtmaConfig;

nan_idx = any(isnan(Measurement),2)' | any(isnan(Kinematics),2)';
calib_idx = TaskStateMasks.use_for_calib>.9;
reject_idx = [0 double(any(diff(TaskStateMasks.target'),2))']; %times when the target changed
arm_not_moving_idx = all(TaskStateMasks.active_assist_weight==0) & all(TaskStateMasks.brain_control_weight==0); %if end-effector was not allowed to move (ie Presentation state) don't use data
result_code = InterpResultCode(TaskStateMasks.result_code);
success_idx = result_code.SuccessfulTrial==rtma.defines.SuccessfulTrial;

% reject for reaction time:
reject_time = 0;	% seconds
fs = 1/rtma.defines.SAMPLE_LENGTH;  % Hz
reject_samples=round(reject_time*fs);
idxs = find(reject_idx>=1);
reject_idx = zeros(size(reject_idx));
if reject_time>0
    for a = 1:length(idxs)
        reject_idx(idxs(a):idxs(a)+reject_samples)=1;
    end
end
% fprintf('Rejecting %d samples at the begining of each state\n',reject_samples)

% square-root transform:
if ~obj.disableSqrt	% SqrtTransform
    Measurement = sqrt(Measurement);
end

% if FilterSigs
% %% Filters:
% *interpolate then apply* or apply separately at each NAN?
% Spike_filter = [0.0932 0.0885 0.0841 0.0799 0.0759 0.0721 0.0685 0.0651 0.0618 0.0587 0.0558 0.0530 0.0503 0.0478 0.0454]'; %from NHP lab - LPF with -3dB cutoff ~1Hz, weighted to more recent values
% Spike_filter = [.0636 .0615 .0594 .0574 .0555 .0536 .0518 .0501 .0484 .0467 .0452 .0436 .0422 .0408 .0394 .0381 .0368 .0355 .0343 .0331 .0321 .0310]'; %JW: new values for 50 Hz sampling
% Kinem_filter = [0.0588 0.0588 0.0588 0.0588 0.0588 0.0588 0.0588 0.0588 0.0588 0.0588 0.0588 0.0588 0.0588 0.0588 0.0588 0.0588 0.0588]'; %from NHP lab - flat LPF window with -3dB cutoff ~1Hz
% Kinem_filter = (1/25).*ones(25,1); % for 50 Hz sampling

t = 1:size(Measurement,1);
if ~isempty (obj.SpikeFilter)
    M =  interp1(t(~nan_idx),Measurement(~nan_idx,:),t,'linear',0);
    % window_len = 0.5*fs;
    % window = hanning(round(window_len));
    % window = window/sum(window);
    % M = conv2(M,window,'same');
    M = conv2(M,obj.SpikeFilter,'valid');
    M = [nan(size(obj.SpikeFilter,1)-1,size(M,2)); M];
    Measurement(~nan_idx,:) = M(~nan_idx,:);
end
if ~isempty(obj.KinematicFilter)
    K =  interp1(t(~nan_idx),Kinematics(~nan_idx,:),t,'linear',0);
    % window_len = 0.5*fs;
    % window = hanning(round(window_len));
    % window = window/sum(window);
    % K = conv2(K,window,'same');
    K = conv2(K,obj.KinematicFilter,'same');
    % K = [nan(size(obj.KinematicFilter,1)-1,size(K,2)); K];
    Kinematics(~nan_idx,:) = K(~nan_idx,:);
end
nan_idx = any(isnan(Measurement),2)' | any(isnan(Kinematics),2)';
% end

% pare down to trainable data:
if ~success_only
    idx = nan_idx | ~calib_idx | reject_idx | arm_not_moving_idx;
else
    idx = nan_idx | ~calib_idx | reject_idx | arm_not_moving_idx | ~success_idx;
end
if remove_bad_trials
    if isfield(TaskStateMasks,'bad_samples')
        idx = idx | TaskStateMasks.bad_samples;
    else
        disp('Cannot reject bad samples, none listed in data file')
    end
end
Measurement(idx,:)=[];
Kinematics(idx,:)=[];

if size(Measurement,1)~=size(Kinematics,1)
    if size(Measurement,1)==size(Kinematics,2)
        Measurement= Measurement';
    else
        weights = [];
        fprintf(2,'Vectors incorrect sizes: Neural(%d , %d), Kinematic(%d , %d) \n',size(Measurement,1),size(Measurement,2),size(Kinematics,1),size(Kinematics,2))
    end
end
num_observations = size(Measurement,1);

% need to make sure all training situations are balanced, if possible  (if
% there are discrete kinematic outputs)
% Kinematics(abs(Kinematics)<0.01*max(max(abs(Kinematics)))) = 0; %sparsify by threshold
utrain = unique(Kinematics,'rows');
if length(utrain)<20
    % find number of each:
    counts = zeros(size(utrain,1),1)';
    for a = 1:size(Kinematics,1)
        which_one(a,:) = all(ones(size(utrain,1),1)*Kinematics(a,:)==utrain,2);
        counts  = counts + which_one(a,:);
    end
    
    max_counts = max(counts);
    for a = 1:size(utrain,1)
        if counts(a)<max_counts
            % randomly pick max_counts-counts(a) observations to repeat
            idx = find(which_one(:,a));
            idxs = unidrnd(length(idx),max_counts-counts(a),1);
            % add onto the training set AND the full set
            Kinematics((end+1):(end+(max_counts-counts(a))),:) = Kinematics(idx(idxs),:);
            Measurement((end+1):(end+(max_counts-counts(a))),:) = Measurement(idx(idxs),:);
        end
    end
end

% standardize / z-score
featureMeans = nanmean(Measurement,1);
featureStds = nanstd(Measurement,1);
obj.featureMeans = featureMeans;
obj.featureStds = featureStds;
% featureStds = ones(size(featureStds));

Measurement_norm = Measurement -  ones(size(Measurement,1),1)*featureMeans; %just for feature selection.
Measurement_norm = Measurement_norm./(ones(size(Measurement_norm,1),1)*featureStds);

%%% FEATURE SELECTION:
if ~isempty(IgnoreKinematicDims) %#ok<NODEF>
    IgnoreKinematicDims(size(Kinematics,2)+1:end) = []; %limit to same size as Kinematic sig, since we are probably not using the last 8 dims.
    Kinematics_masked  = Kinematics(:,IgnoreKinematicDims==0);
else
    Kinematics_masked = Kinematics;
end

obj.MeasurementMask = obj.FeatureSelection(Measurement_norm,Kinematics_masked,0.005,0,obj.PreTrainMask);

end

% mask(mean(Measurement,1)<.01)=0; % minimum spike rate
% mask(1:5:end) = 0; % no hash
