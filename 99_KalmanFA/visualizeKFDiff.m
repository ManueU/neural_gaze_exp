% This is a function to visualize the parameters of two different kalman
% filters. 
function visualizeKFDiff(kf1, kf2)

f = figure(); 
f.Units = 'normalized';

% Compare state models 
statePanel = uipanel('Parent', f, 'Position', [0 .5 1 .5]); 

mu1Ax = axes('Parent', statePanel, 'Position', [.05 .1 .15 .7]); 
plot(mu1Ax, kf1.mu_1, 'bo');
hold on; 
plot(mu1Ax, kf2.mu_1, 'r*'); 
hold off;
title('Parent', mu1Ax, 'mu_1');

v1Panel = uipanel('Parent', statePanel, 'Position', [.25 0 .25 1]); 
plotNMats({kf1.V_1, kf2.V_1, kf2.V_1 - kf1.V_1}, {'V1: K1', 'V1: K2', 'V1: K2 - K1'}, 'LAYOUT', [3, 1], 'PARENT', v1Panel, 'FONT_SIZE', 10); 

aPanel = uipanel('Parent', statePanel, 'Position', [.5 0 .25 1]); 
plotNMats({kf1.A, kf2.A, kf2.A - kf1.A}, {'A: K1', 'A: K2', 'A: K2 - K1'}, 'LAYOUT', [3, 1], 'PARENT', aPanel, 'FONT_SIZE', 10);

qPanel = uipanel('Parent', statePanel, 'Position', [.75 0 .25 1]); 
plotNMats({kf1.Q, kf2.Q, kf2.Q - kf1.Q}, {'Q: K1', 'Q: K2', 'Q: K2 - K1'}, 'LAYOUT', [3, 1], 'PARENT', qPanel, 'FONT_SIZE', 10);

% Compare observation models
obsPanel = uipanel('Parent', f, 'Position', [0 0 1 .5]); 

dAx = axes('Parent', obsPanel, 'Position', [.05 .1 .15 .7]); 
plot(dAx, kf1.d, 'bo');
hold on; 
plot(dAx, kf2.d, 'r*'); 
hold off;
title('Parent', dAx, 'd');

cPanel = uipanel('Parent', obsPanel, 'Position', [.25 0 .5 1]); 
plotNMats({kf1.C, kf2.C, kf2.C - kf1.C}, {'C: K1', 'C: K2', 'C: K2 - K1'}, 'LAYOUT', [1, 3], 'PARENT', cPanel, 'FONT_SIZE', 10);
colorbar(); 

rPanel = uipanel('Parent', obsPanel, 'Position', [.75 0 .25 1]); 
plotNMats({kf1.R, kf2.R, kf2.R - kf1.R}, {'R: K1', 'R: K2', 'R: K2 - K1'}, 'LAYOUT', [1, 3], 'PARENT', rPanel, 'FONT_SIZE', 10);
colorbar(); 


