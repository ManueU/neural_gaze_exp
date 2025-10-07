% This is a script to test our Kalman Filter code

%% ========================================================================
%                        PARAMETERS GO HERE
%  ========================================================================
nSteps = poissrnd(5, [1000,1]);

xDim = 2;
yDim = 10;

trueKf.d = randn(yDim, 1)*10;
trueKf.A = eye(xDim);
trueKf.C = randn(yDim, xDim);
trueKf.Q = eye(xDim);
trueKf.R = generateRandomCovMatrix(yDim);
trueKf.mu_1 = randn(xDim,1)*100;
trueKf.V_1 = 1000*generateRandomCovMatrix(xDim);

%% ========================================================================
%                  GENERATE RANDOM TRIALS FROM THE TRUE MODEL
%  ========================================================================
[trueX, trueY] = generateKalmanData(trueKf, nSteps);

%% ========================================================================
%                   DECODE TRIALS WITH THE TRUE MODEL
%  ========================================================================
[xHat, postCov] = kalmanDecode(trueKf, trueY);

%% ========================================================================
%                           VISUALIZE A TRIAL
%  ========================================================================

figure();
plotAx = axes();
for pT = 1:length(nSteps)
    if ~isempty(trueX{pT})
        plot(trueX{pT}(1,:), trueX{pT}(2,:), 'b--', 'lineWidth', 2);
        hold on;
        plot(xHat{pT}(1,:), xHat{pT}(2,:), 'r-', 'lineWidth', 2);
        
        curNSteps = size(xHat{pT},2);
        for s = 1:1:curNSteps
            pH = plot2DNormalConfEllipse(xHat{pT}(:,s), postCov{pT}(:,:,s), .95, 100, 'AX', plotAx);
            pH.Color = [1 0 0];
            hold on;
        end
        
    end
    hold off;
    pause
end
%% ========================================================================
%                      FIT THE KALMAN FILTER TO DATA
%  ========================================================================
estKf = kalmanFit(trueX,trueY, 'CALC_CONV_GAIN', true);

%% ========================================================================
%                   VISUALIZE THE FIT KALMAN FILTERS
%  ========================================================================
visualizeKFDiff(trueKf, estKf)