% We consider a Kalman filter (KF) for the following dynamical system:
%
%   x_{t+1} = Ax_t + w_t
%       y_t = Cx_t + d + v_t ,
%
% where x_t is a latent state and y_t is an observed state, and we make the
% following distributional assumptions:
%
%   x_0 ~ N(mu_1, V_1)
%   w_t ~ N(0, Q)
%   y_t ~ N(0, R)
%
% Given a structure with the parameters of this model, this function will
% generate random trials according to the model.
%
% Usage: [x, y]  = generateKalmanData(kf, nSteps, varargin)
%
% Input:
%
%   kf - a structure with the parameters of the dynamical model for the
%   filter.  It should contain the fields d, A, C, mu_1, V_1, Q and R.
%
%   nSteps - an array length N, where N is the number of trials to simulate.
%          nSteps(n) should contain the number of steps to simulate for
%          trial n.
%
% Optional Inputs: All optional inputs should be entered in string-value
% pair format. 
%
%   NO_ASSERT - Assertions are skipped.  Good for efficiency but may not
%               catch unexpected inputs. Default: true
%           
%   VERBOSE - True if progress updates should be output to screen.  
%             Default: true
%
% Output:
%
%       x - a cell of length N.  Each entry contains a the latent state for a trial as a 
%	matrix of size xDim by T, where xDim is the dimensionality of x and T is the
%	number of steps in a particular trial
%
%       y - a cell of length N.  Each entry contains the observations for the corresponding 
%	trial in x as a  matrix of size yDim by T, where yDim is the dimensionality of y.
%
% Author: wbishop@cs.cmu.edu
%
function [x,y] = generateKalmanData(kf, nSteps, varargin)

NO_ASSERT = false;
VERBOSE = true; 
warnOpts(assignOpts(varargin)); 

if ~NO_ASSERT
    verifyKF(kf); 
    assert(all(nSteps >= 0), 'Found negative step counts');
end

nTrials = length(nSteps);
xDim = length(kf.A); 
yDim = size(kf.C,1); 

% Generate initial states for all trials
x0 = mvnrnd(kf.mu_1', kf.V_1, nTrials)'; 

% Simulate each trial
x  = cell(1, nTrials); 
y = cell(1, nTrials); 
for t = 1:nTrials
    curNSteps = nSteps(t); 
    if curNSteps > 0
        stateTransitionNoise = mvnrnd(zeros(1, xDim), kf.Q, curNSteps-1)';
        obsNoise = mvnrnd(zeros(1, yDim), kf.R, curNSteps)';
        
        curX = nan(xDim, curNSteps);
        curY = nan(yDim, curNSteps);
        curX(:,1) = x0(:,t);
        curY(:,1) = kf.C*curX(:,1) + kf.d + obsNoise(:,1);
        for s = 2:curNSteps
            curX(:,s) = kf.A*curX(:,s-1) + stateTransitionNoise(:,s-1);
            curY(:,s) = kf.C*curX(:,s) + kf.d + obsNoise(:,s);
        end
    else
        curX = []; 
        curY = []; 
    end
    x{t} = curX;  
    y{t} = curY;  
    
    if VERBOSE && mod(t,100) == 0
        disp(['Completed simulating trial ', num2str(t) ' of ', num2str(nTrials), '.']); 
    end
end

