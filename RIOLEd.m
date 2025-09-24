
classdef RIOLEd < Decoders.Decoder
    %Decoders.RIOLEd
    % Ridge Regression Indirect Optimal Linear Estimator decoder
    
    properties (SetAccess = immutable)
       version = 1.0; 
    end
    
    properties
       Kmeans % consider moving to Decoders.Decoder
       Kstd   % consider moving to Decoders.Decoder
       baselines
       weights
    end
    
    methods
        function obj = RIOLEd(varargin)
            if nargin > 0
               obj = obj.SetConfig(varargin{:}); 
            end
        end
        
    end
    
end

