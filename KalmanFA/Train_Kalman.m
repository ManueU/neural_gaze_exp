function [obj,M2] = Train(obj,data_files,opt)
global RTMA

% [D,M2] = Train(D,NeuralSig,KinematicSig,NeuralSig_val,KinematicSig_val)
%
% Train method for the KalmanFA decoder class.
%
% Decoding parameters for the KalmanFA decoder are computed offline.  This
% function serves as a placeholder that can be used to load a set of
% previously-computed decoding parameters into memory.
%
% Normally this function will not be called.  Instead, a decoder object
% will be trained offline and loaded from disk directly.
%
% Inputs:
%   D       KalmanFA decoder object
%   y       Spike count vector
%
% Outputs:
%   D       KalmanFA decoder object
%   M2      M2 weights matrix.
%
% Author:       Alan D. Degenhart
% Date Created: 2016.05.19
% Last Updated: 2016.05.19
% Load decoder

[f,p] = uigetfile('*.mat','Select Decoder Weights');
K = load([p f]);
K = K.kalmanInitParams;

% Put weights into decoder object
obj.M0 = K.M0;
obj.M1 = K.M1;
obj.M2 = K.M2;
obj.mu0 = K.mu0;

end

