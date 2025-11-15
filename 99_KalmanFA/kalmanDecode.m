% A function to calculate P(x_t| {y}_1^t) using a Kalman filter
%
% Usage: [xHat, postCov] = kalmanDecode(kf, data, varargin)
%
% Inputs:
%
%   kf - a structure of parameters for the Kalman filter.  See
%   generateKalmanData for more information.
%
%   data - a cell of length N, where N is the number of trials.  data{n}
%   contains observed data, Y, for the n^th trial.  Y should be size yDim by T,
%   where yDim is the observation dimensionality and T is the number of steps
%   in the trial.
%
% Optional Inputs: All optional inputs should be given in string-value pair
% format.
%
%   NO_ASSERT - Assertions are skipped.  Good for efficiency but may not
%               catch unexpected inputs. Default: false
%
%   VERBOSE - True if progress updates should be output to screen.
%             Default: true
%
% Outputs:
%
%       xHat - a cell of length N.  Each entry contains a matrix of size
%	xDim by T, where xDim is the dimensionality of the latent state for 
%	the Kalman filter.  xHat{n}(:,t) contains the posterior mean 
%	for P(x_t| {y}_t^1) for the n^th trial.
%
%       postCov - (Optional) A cell of length N.  Each entry contains an array
%	of size xDim  by XDim  by T.  postCov{n}(:,:,t) contains the covariance for 
%	P(x_t| {y}_t^1) for the n^th trial. 
%
% Author: wbishop@cs.cmu.edu
%

function [xHat, postCov]  = kalmanDecode(kf, data, varargin)

NO_ASSERT = false;
VERBOSE = true;
warnOpts(assignOpts(varargin));

if ~NO_ASSERT
    verifyKF(kf);
end

recordCov = nargout > 1; 

nTrials = length(data);
xDim = length(kf.A);


V_1 = kf.V_1;
mu_1 = kf.mu_1;
A = kf.A;
C = kf.C;
R = kf.R;
Q = kf.Q;
d = kf.d;

xHat = cell(1, nTrials); 
if recordCov
    postCov = cell(1, nTrials); 
end
for t = 1:nTrials
    curY = data{t};
    if ~isempty(curY)
        curYCtr = bsxfun(@minus, curY, d);
        curNSteps = size(curY,2);
        
        curXHat = nan(xDim, curNSteps);
        if recordCov
            curPostCov = nan(xDim, xDim, curNSteps);
        end
        
        for s = 1:curNSteps
            if s == 1
                oneStepV = V_1;
                oneStepX = mu_1;
            else
                oneStepV = A*prevV*A' + Q;
                oneStepX = A*curXHat(:,s-1);
            end
            
            curK = oneStepV*C'*inv(R + C*oneStepV*C');
            curV = oneStepV - curK*C*oneStepV;
            curXHat(:, s) = oneStepX + curK*(curYCtr(:,s) - C*oneStepX);
            
            if recordCov
                curPostCov(:, :, s) = curV;
            end
            prevV = curV;
        end
    else
        curXHat = []; 
        curPostCov = []; 
    end
    
    xHat{t}= curXHat;
    if recordCov
        postCov{t} = curPostCov;
    end
    
    if VERBOSE && mod(t,100) == 0;
        disp(['Done decoding trial ', num2str(t), ' of ', num2str(nTrials)]);
    end
end
