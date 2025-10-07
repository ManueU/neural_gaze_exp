clear 
% close all
clc


load('controlled_BCI02.mat')
n_sets = 6; 
n_trials = 32*ones(1,n_sets);  
bin_size = 0.02;
period_pre = 0.1; 
period_reach = 0.5; 


%% Data extraction
% Target IDs
tmp = []; 
for set = 1:n_sets
    Y = [tmp; [data(set).Data(1).Resampled.Target_ID]']; 
    tmp = Y;
end 


% Firing rate  
idx_pres = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == "Pres12"); 
idx_reach = find(string(data(1).Data(1).Resampled(1).Task_states(:,1)) == "Reach"); 

start_pres = size(data(1).Data(2).Resampled(1).Task_states{idx_pres,2}, 1) - (period_pre/bin_size); 
end_reach = period_reach/bin_size;

j = 1; 
X = cell(sum(n_trials),1);
for set = 1:n_sets
    for trial = 1:n_trials(set)
        tmp_pres = []; 
        tmp_reach = []; 
        for array = 1:2
            tmp_pres = [tmp_pres, data(set).Data(array).Resampled(trial).Task_states{idx_pres,2}(start_pres:end, :)]; 
            tmp_reach = [tmp_reach, data(set).Data(array).Resampled(trial).Task_states{idx_reach,2}(1:end_reach, :)]; 
        end 
        matrix = [tmp_pres; tmp_reach]; 
        X{j} = mean(matrix./bin_size,1);
        j = j + 1; 
    end   
end 
X = cell2mat(X); 

%% Divisione training/test
cv = cvpartition(Y,'HoldOut',0.3,'Stratify',true);
idxTrain = training(cv); 
idxTest = test(cv);

Xtrain = X(idxTrain, :);   
Xtest = X(idxTest, :);   

Ytrain = Y(idxTrain);
Ytest  = Y(idxTest);


%% Addestramento decoder (SVM multiclass)
t = templateSVM('KernelFunction','rbf','KernelScale','auto','Standardize',true);
Msvm = fitcecoc( ...
    Xtrain, Ytrain, ...
    'Learners', t, ...
    'Coding', 'onevsall');

%% Predizione
Ypred = predict(Msvm, Xtest);

% Accuratezza
accuracy = mean(Ypred == Ytest);
disp(['Decoder accuracy: ', num2str(accuracy*100, '%.2f'), '%']);

% Confusion matrix con percentuali
classes = unique(Y);
classLabels = arrayfun(@(c) sprintf('Target %d', c), classes, 'UniformOutput', false);
Ytest_cat = categorical(Ytest, 1:length(classLabels), classLabels);
Ypred_cat = categorical(Ypred, 1:length(classLabels), classLabels);
figure('Color','w');
cm = confusionchart(Ytest_cat, Ypred_cat);
cm.Normalization = 'row-normalized';  % percentuali su ogni riga
cm.Title = sprintf('Confusion Matrix - Location Decoder (Accuracy: %.2f%%)', accuracy*100);
cm.XLabel = 'Predicted Target';
cm.YLabel = 'True Target';

