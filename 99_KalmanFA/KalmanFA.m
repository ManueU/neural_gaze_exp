classdef KalmanFA < Decoders.Decoder
    % Author:       Alan D. Degenhart
    % Date Created: 2016.05.19
    % Last Updated: 2018.07.11
    
    properties (SetAccess = immutable)
        version = 1.0;
    end
    
    properties
        M0      % Offset vector
        M1      % Influence of previous state estimate
        M2      % Influence of neural activity
        mu0     % Initial state estimate
        x_prev	% Previous state estimate
        weights	% Weights vector (visualization ONLY)
        y_sum	% Previous spike count (to make up 80 ms bin)
        counter	% Counter to increase bin size
    end
    
    methods
        function obj = KalmanFA(Config)
            if exist('Cofig','var')
                obj.SetConfig(Config);
            end
        end
    end
end

