%% MLP DECODER — Neural / Gaze / Neural+Gaze

fprintf('\n================ MLP DECODER ================\n');

%% Standardizzazione input/output
[X_train_z, X_val_z, muX, sigX] = standardize(X_train_decoder, X_val_decoder);
X_test_z = (X_test_decoder - muX) ./ sigX;

[Y_train_z, Y_val_z, muY, sigY] = standardize(Y_train, Y_val);

X_train_z(isnan(X_train_z) | isinf(X_train_z)) = 0;
X_val_z(isnan(X_val_z) | isinf(X_val_z)) = 0;
X_test_z(isnan(X_test_z) | isinf(X_test_z)) = 0;

%% MATLAB neural network vuole: features x samples
Xtr = X_train_z';
Ytr = Y_train_z';

Xva = X_val_z';
Yva = Y_val_z';

Xte = X_test_z';

%% Architettura MLP piccola
hidden_units = 16;

net = fitnet(hidden_units, 'trainscg');  % scaled conjugate gradient

net.performFcn = 'mse';

net.layers{1}.transferFcn = 'poslin';   % ReLU
net.layers{2}.transferFcn = 'purelin';  % output lineare

%% Divisione train/validation manuale
net.divideFcn = 'divideind';

n_train = size(Xtr, 2);
n_val   = size(Xva, 2);

X_all = [Xtr, Xva];
Y_all = [Ytr, Yva];

net.divideParam.trainInd = 1:n_train;
net.divideParam.valInd   = n_train + (1:n_val);
net.divideParam.testInd  = [];

%% Regolarizzazione / early stopping
net.trainParam.max_fail = 10;
net.trainParam.epochs   = 500;
net.trainParam.showWindow = true;

net.performParam.regularization = 0.2;

%% Train
[net, tr] = train(net, X_all, Y_all);

%% Predict validation
Y_pred_val_z = net(Xva)';
Y_pred_val = Y_pred_val_z .* sigY + muY;
Y_val_eval = Y_val;

RMSE_val = sqrt(mean((Y_val_eval - Y_pred_val).^2, 1));

SS_res_val = sum((Y_val_eval - Y_pred_val).^2, 'all');
SS_tot_val = sum((Y_val_eval - mean(Y_val_eval, 1)).^2, 'all');
R2_val = 1 - SS_res_val / SS_tot_val;

fprintf('\n--- Validation performance MLP ---\n');
fprintf('R² val globale: %.4f\n', R2_val);
fprintf('RMSE vx val: %.4f | vy val: %.4f\n', RMSE_val(1), RMSE_val(2));

%% Predict test
Y_pred_test_z = net(Xte)';
Y_pred_test = Y_pred_test_z .* sigY + muY;
Y_test_eval = Y_test;

%% Metriche test
SS_res_test = sum((Y_test_eval - Y_pred_test).^2, 'all');
SS_tot_test = sum((Y_test_eval - mean(Y_test_eval, 1)).^2, 'all');

R2_test = 1 - SS_res_test / SS_tot_test;

R2_test_dim = 1 - ...
    sum((Y_test_eval - Y_pred_test).^2, 1) ./ ...
    sum((Y_test_eval - mean(Y_test_eval, 1)).^2, 1);

RMSE_test = sqrt(mean((Y_test_eval - Y_pred_test).^2, 1));

fprintf('\n--- Test performance MLP ---\n');
fprintf('R² globale TEST : %.4f\n', R2_test);
fprintf('R² vx TEST      : %.4f\n', R2_test_dim(1));
fprintf('R² vy TEST      : %.4f\n', R2_test_dim(2));
fprintf('RMSE vx TEST    : %.4f\n', RMSE_test(1));
fprintf('RMSE vy TEST    : %.4f\n', RMSE_test(2));

%% Plot test
figure('Color','w');

subplot(2,2,1);
plot(Y_test_eval(:,1), 'k'); hold on;
plot(Y_pred_test(:,1), 'r--');
xlabel('Time bin'); ylabel('v_x');
title(sprintf('MLP TEST vx — R²=%.3f', R2_test_dim(1)));
legend('Reale','Predetta'); grid on;

subplot(2,2,2);
plot(Y_test_eval(:,2), 'k'); hold on;
plot(Y_pred_test(:,2), 'r--');
xlabel('Time bin'); ylabel('v_y');
title(sprintf('MLP TEST vy — R²=%.3f', R2_test_dim(2)));
legend('Reale','Predetta'); grid on;

subplot(2,2,3);
scatter(Y_test_eval(:,1), Y_pred_test(:,1), 8, 'filled'); hold on;
refline(1,0);
xlabel('vx reale'); ylabel('vx predetta');
title('MLP scatter vx'); grid on;

subplot(2,2,4);
scatter(Y_test_eval(:,2), Y_pred_test(:,2), 8, 'filled'); hold on;
refline(1,0);
xlabel('vy reale'); ylabel('vy predetta');
title('MLP scatter vy'); grid on;

%% Plot vx-vy trial-by-trial
figure('Color','w');
hold on;

test_trial_ids = unique(trial_id_test, 'stable');

for i = 1:numel(test_trial_ids)

    tr_id = test_trial_ids(i);
    idx = trial_id_test == tr_id;

    plot(Y_test_eval(idx,1), Y_test_eval(idx,2), ...
        'k-', 'LineWidth', 1.5);

    plot(Y_pred_test(idx,1), Y_pred_test(idx,2), ...
        'r--', 'LineWidth', 1.5);
end

xlabel('v_x');
ylabel('v_y');
title('MLP TEST — vx-vy trial-by-trial');

axis equal;
grid on;
legend('Reale','Predetta');