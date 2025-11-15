function [obj, ProcessStateEst] = Predict( obj, Measurement)

% [OLE, ProcessStateEst] = Predict( OLE, Measurement)
%
% Perform an update/predict step for the OLE Decoder. Measurement
% is the observed data. The output ProcessStateEst represents the estimated
% state of the process. No smoothing is performed on the Measurement. If
% smoothing is desired, it should be done externally prior to passing the
% Measurement to this function.

% Brian Wodlinger based on Meel Velliste 9/19/2009

Measurement = Measurement -  ones(size(Measurement,1),1)*obj.featureMeans;
Measurement = Measurement./(ones(size(Measurement,1),1)*obj.featureStds);


% Grab stuff from config
if ~isempty(obj.MeasurementMask)
	MMask = logical(obj.MeasurementMask);
else
	MMask = true(1,size(Measurement,2));
	%fprintf('No Mask available\n')
end
if isempty(obj.weights)
	ProcessStateEst = nan(1,size(obj.weights,2));
	fprintf('No weights available\n')
	return;
end
weights = obj.weights;

Measurement = Measurement(:,MMask);

if size(Measurement,2)~=size(weights,1)
	if size(Measurement,1)==size(weights,1)
		Measurement = Measurement';
	else
		ProcessStateEst = nan(1,size(weights,2));
		fprintf(2,'Measurement vector incorrect size: %d %d (expected %d)\n',size(Measurement,1),size(Measurement,2),size(weights,1)-1)
        return
	end
end
 
% Estimate current process state based on measurement
ProcessStateEst = (Measurement-repmat(obj.baselines,[size(Measurement,1) 1]))*weights;

if ~isempty(obj.Kmeans)
    ProcessStateEst = ProcessStateEst.*(ones(size(ProcessStateEst,1),1)*obj.Kstd);
    ProcessStateEst = ProcessStateEst + (ones(size(ProcessStateEst,1),1)*obj.Kmeans);
end

end
