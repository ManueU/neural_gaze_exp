% This is a function to run a set of assertions on a Kalman filter
% structure. 
%
% Usage: verifyKF(kf)
%
% Inputs: 
%
%   kf - a Kalman filter structure (see function generateKalmanData)
%
% Outputs: 
%
%   None. 
%
% Author: wbishop@cs.cmu.edu
%
function verifyKF(kf)

xDim = length(kf.A); 
yDim = size(kf.C,1); 

assert(size(kf.A,1) == size(kf.A,2), 'A must be a square matrix.');

assert(ismatrix(kf.d), 'd must be a column vector.'); 
assert(size(kf.d,2) == 1, 'd must be a column vector.'); 
assert(size(kf.d,1) == yDim, 'd is the wrong sixe for C');

assert(size(kf.C,2) == xDim, 'C is the wrong size for A.'); 

assert(size(kf.mu_1,2) == 1 && ismatrix(kf.mu_1), 'mu_1 must be a column vector'); 
assert(size(kf.mu_1,1) == xDim, 'mu_1 is the wrong size for A.'); 

assert(ismatrix(kf.V_1), 'V_1 must be a covariance matrix.'); 
assert(size(kf.V_1,1) == size(kf.V_1,2), 'V_1 must be a covariance matrix.'); 

assert(ismatrix(kf.Q), 'Q must be a covariance matrix.'); 
assert(size(kf.Q,1) == size(kf.Q,2), 'Q must be a covariance matrix.'); 
assert(size(kf.Q,1) == xDim, 'Q is the wrong size for A'); 

assert(ismatrix(kf.R), 'R must be a covariance matrix.'); 
assert(size(kf.R,1) == size(kf.R,2), 'R must be a covariance matrix.'); 
assert(size(kf.R,1) == yDim, 'R is the wrong size for C'); 
