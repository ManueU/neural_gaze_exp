clear all
close all
clc

[file, path] = uigetfile('*.mat', 'Seleziona un file MAT');
if isequal(file,0)
    error('Nessun file selezionato.');
else
    disp('Loading data...')
    load(fullfile(path, file));
end

%% Input data
% Spike count matrix
MaskedSpikeCounts = set01.SpikeCount(:,1:5:end);
Motor1SC = MaskedSpikeCounts(:,[1:64 97:128]);
Stim1SC = MaskedSpikeCounts(:,65:96);
Motor2SC = MaskedSpikeCounts(:,[129:192 225:256]);
Stim2SC = MaskedSpikeCounts(:,193:224);
motor_data = [Motor1SC, Motor2SC]; 

% Velocities in output from the neural decoder
vel(:,1) = set01.Kinematics.Control(:,2); 
vel(:,2) = set01.Kinematics.Control(:,3); 
pos(:,1) = set01.Kinematics.Actual(:,2); 
pos(:,2) = set01.Kinematics.Actual(:,3); 
cursor_kin = [pos, vel];

% Target_ID labeling 
target_coordinates = [set01.TaskStateMasks.target(2, :)', set01.TaskStateMasks.target(3, :)'];
target_info = target_coordinates; 
target_info(:,3) = NaN;
refXY = [  0      0
           0      0.2
           0.14   0.14
           0.2    0
           0.14  -0.14
           0     -0.2
          -0.14  -0.14
          -0.2    0
          -0.14   0.14 ];
codes = (0:8).';
tol = 1e-6;
[tf, loc] = ismembertol(target_info(:,1:2), refXY, tol, 'ByRows', true);
target_info(tf,3) = codes(loc(tf));
target_id = target_info(:,3);
target_id(isnan(target_id)) = 0; 

% Phase labeling
phase = set01.TaskStateMasks.state_name; 
phase(cellfun('isempty', phase)) = {'Center'}; 

 
