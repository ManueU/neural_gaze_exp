function EMConfig = convertKalmantoEM(K,DATALOC)
% [D] = convertKalmanParams(K)
%
% Conver Kalman filter params to HST KalmanFA decoder object.
%
% This function takes the parameters in the 'kalmanInitParams' Kalman
% filter decoder structure and converts them to a 'KalmanFA' decoder object
% compatible with the HST system.
%
% Inputs:
%   K       'kalmanInitParams' decoding parameter structure
%
% Outputs:
%   D       'KalmanFA' decoder object
%
% Author:       Alan D. Degenhart
% Date Created: 2016.05.20
% Last Updated: 2016.05.20

% Get appropriate parameters

C.M0 = K.M0;                         % Offset vector
C.M1 = K.M1;                         % Influence of previous state estimate
C.M2 = K.M2;                         % Influence of neural activity
C.mu0 = K.mu0;                       % Initial state estimate
C.x_prev = zeros(length(DATALOC),1); % Previous state estimate

% Get feature mask and initialize feature mean and standard deviation
% vectors
chMask = K.validChMask;
nUnits = length(chMask);
featMean = zeros(1,nUnits);
featStd = ones(1,nUnits);

% Get feature mean and standard deviation from kalmanInitParams object
if ~isempty(K.NormalizeSpikes.mean)
    featMean(chMask) = K.NormalizeSpikes.mean;
    featStd(chMask) = K.NormalizeSpikes.std;
end

% Set weights
nUnits = length(chMask); %changed from sum(chMask) because chMask no longer ones vector
weights = zeros(nUnits,30);
weights(:,DATALOC) = K.M2(:,chMask)';
C.weights = weights;  % Weights vector (visualization ONLY)

% Update KalmanFA object with parameters
C.MeasurementMask = chMask;      % Indicates which features are used
C.featureMeans = featMean;       % Mean of features
C.featureStds = featStd;         % Standard deviation of features
C.y_sum = zeros(1,1280);         % Previous spike count (to make up 80 ms bin)
C.counter = 1;                   % Counter to increase bin size

% Instantiate KalmanFA decoder object with specified parameters
D = KalmanFA(C);

% Convert to EMConfig object
EMConfig.DecoderType = 'KalmanFA';
EMConfig.DecoderConfig = [];
EMConfig.Decoder = D;
% EMConfig.pertParams = K.P.ps.pertParams;
if isfield(K.P,'psStitch'); EMConfig.decoderParams = K.P.psStitch; end
EMConfig.psBase = K.P.psBase;
EMConfig.decoderNum = K.P.dailyDecoderNum;
