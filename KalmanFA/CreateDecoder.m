function Decoder = CreateDecoder(DecoderType,DecoderConfig)
%
% Create decoder object and set the initial config
%

if ~exist('DecoderConfig','var')
    DecoderConfig = [];
end

% create function handle (assuming full name was input)
hDecoder = str2func(DecoderType);
try
    Decoder = hDecoder(DecoderConfig);
catch
    % try a short name instead
    switch( DecoderType)
        case {'pva','PopulationVector'}
            Decoder = PopulationVector( DecoderConfig);
            
        case {'kf','KalmanFilter'}
            Decoder = KalmanFilter( DecoderConfig);
            
        case {'lgf','LaplaceGaussianFilter'}
            Decoder = LaplaceGaussianFilter( DecoderConfig);
            
        case {'ole','OptimalLinearEst'}
            Decoder = OptimalLinearEst( DecoderConfig);
            
        case {'iole','InvOptimalLinearEst'}
            Decoder = InvOptimalLinearEst( DecoderConfig);
            
        case {'svm','SupportVectorMachine'}
            Decoder = SupportVectorMachine( DecoderConfig);
            
        case {'nhpole','NHPOLE'}
            Decoder = NHPOLE( DecoderConfig);
            
        case {'role','ROptimalLinearEst'}
            Decoder = ROptimalLinearEst( DecoderConfig);
            
        case {'riole','RIndOptimalLinearEst'}
            Decoder = RIndOptimalLinearEst( DecoderConfig);
        case {'riole1','RIndOptimalLinearEst1'}
            Decoder = RIndOptimalLinearEst1( DecoderConfig);
            
        case {'aole','AdaptiveOLE'}
            Decoder = AdaptiveOLE( DecoderConfig);
            
        case {'rf','RandomForest'}
            Decoder = RandomForest( DecoderConfig);

        case {'RIOLEd'}
            Decoder = Decoders.RIOLEd( DecoderConfig);
        case {'RIOLEdForceVelocity'}
            Decoder = RIOLEdForceVelocity( DecoderConfig);
        case {'LDA'}
            Decoder = LDA( DecoderConfig);
        otherwise
            error( ['Unrecognized Decoder type ''' DecoderType '''!']);
    end
    
end

end
