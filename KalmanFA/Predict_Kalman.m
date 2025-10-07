function [obj,x] = Predict(obj,y)
% [obj,x] = Predict(obj,y)
%
% Predict method for the KalmanFA decoder class.
%
% Inputs:
%   obj       KalmanFA decoder object
%   y         Spike count vector
%
% Outputs:
%   obj       KalmanFA decoder object
%   x         Decoded state estimate (generally velocity)
%
% Author:       Alan D. Degenhart
% Date Created: 2016.05.19
% Last Updated: 2018.07.20 by AH to update combined bin logic (counter)

% Get decoding parameters
M0 = obj.M0;           % M0 matrix
M1 = obj.M1;           % M1 matrix
M2 = obj.M2;           % M2 matrix
x_prev = obj.x_prev;   % State est. from previous timestep
% y_sum = obj.Members.Config.y_sum;     % Previous spike count (to make up 80 ms bin)
% counter = obj.Members.Config.counter; % Counter to increase bin size

% Add a counter and logic to only update after counter reaches 4 (to get 80
% ms bins)
% Predict state estimate
x = M0 + M1*x_prev + M2*y';
% if counter == 3
%     x = M0 + M1*x_prev + M2*(y_sum + y)';
%     y_sum = zeros(1,1280);
%     counter = 1;
% else
%     x = x_prev;
%     y_sum = y_sum + y;
%     counter = counter+1;
% end

% Put state estimate into decoder object for next time step
obj.x_prev = x;
% obj.Members.Config.y_sum = y_sum;
% obj.Members.Config.counter = counter;

% Zero-pad state estimate to the size expected by the HST system (1 x 30)
x_full = zeros(1,30);
x_full(2:3) = x';
x = x_full;

end

