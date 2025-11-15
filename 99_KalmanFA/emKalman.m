function kf = emKalman(yCell, xDim, maxIters, tol)
    if nargin < 3, maxIters = 100; end
    if nargin < 4, tol = 1e-6; end
    
    % --- Preparazione ---
    N = numel(yCell);
    yDim = size(yCell{1},1);
    Tsum = 0; for n=1:N, Tsum = Tsum + size(yCell{n},2); end
    
    % Inizializzazione semplice
    Yall = []; for n=1:N, Yall = [Yall, yCell{n}]; end
    d = mean(Yall, 2);
    Yc = Yall - d;
    [U,S,~] = svd(Yc, 'econ');  % PCA
    C = U(:,1:xDim)';
    C = C';                     % yDim x xDim
    Xinit = C \ Yc;             % LS sugli scores
    A = eye(xDim);
    if size(Xinit,2) > 1
        Xp = Xinit(:,2:end); Xm = Xinit(:,1:end-1);
        A = Xp * Xm' / (Xm*Xm' + 1e-6*eye(xDim));
    end
    R = cov((Yc - C*Xinit)'); if isempty(R), R = eye(yDim); end
    Q = eye(xDim) * 1e-2;
    mu_1 = zeros(xDim,1);
    V_1 = eye(xDim);
    
    ll_prev = -inf;
    
    for iter = 1:maxIters
        % Accumulatori M-step
        Syy = zeros(yDim); Syx = zeros(yDim,xDim); Sxx = zeros(xDim); % per C,R
        Sx1  = zeros(xDim,1); Sx1x1 = zeros(xDim);
        Sxxtm1 = zeros(xDim); Sxtm1tm1 = zeros(xDim); Sxxt = zeros(xDim); % per A,Q
        Ttot = 0; TtotA = 0;
        ll = 0;
    
        for n = 1:N
            Y = yCell{n};
            Tn = size(Y,2); if Tn==0, continue; end
            Ttot = Ttot + Tn; TtotA = TtotA + max(0,Tn-1);
    
            % --- Kalman forward ---
            x_f = zeros(xDim,Tn); P_f = zeros(xDim,xDim,Tn);
            x_p = zeros(xDim,Tn); P_p = zeros(xDim,xDim,Tn);
            Kt  = zeros(xDim,yDim,Tn); St = zeros(yDim,yDim,Tn);
    
            for t=1:Tn
                if t==1
                    xpred = mu_1; Ppred = V_1;
                else
                    xpred = A * x_f(:,t-1);
                    Ppred = A * P_f(:,:,t-1) * A' + Q;
                end
                yc = Y(:,t) - d - C * xpred;
                S = C*Ppred*C' + R;
                K = Ppred*C'/S;
                x_f(:,t) = xpred + K*yc;
                P_f(:,:,t) = Ppred - K*C*Ppred;
    
                x_p(:,t) = xpred; P_p(:,:,t) = Ppred;
                Kt(:,:,t) = K; St(:,:,t) = S;
                % log-like contrib
                ll = ll - 0.5*( log(det(S)) + yc'/S*yc + yDim*log(2*pi) );
            end
    
            % --- RTS smoother ---
            x_s = x_f; P_s = P_f;
            Jt = zeros(xDim,xDim,Tn-1);
            for t=Tn-1:-1:1
                J = P_f(:,:,t) * A' / P_p(:,:,t+1);
                x_s(:,t) = x_f(:,t) + J*(x_s(:,t+1) - x_p(:,t+1));
                P_s(:,:,t) = P_f(:,:,t) + J*(P_s(:,:,t+1) - P_p(:,:,t+1))*J';
                Jt(:,:,t) = J;
            end
    
            % Cross-covarianze smoothed (forma standard a coeff. costanti)
            P_sm_cross = zeros(xDim,xDim,Tn-1);
            for t=Tn-1:-1:1
                P_sm_cross(:,:,t) = Jt(:,:,t) * P_s(:,:,t+1);
            end
    
            % --- Accumuli per M-step ---
            % Osservazione
            for t=1:Tn
                Syy = Syy + (Y(:,t)-d)*(Y(:,t)-d)';
                Syx = Syx + (Y(:,t)-d)*x_s(:,t)';
                Sxx = Sxx + (P_s(:,:,t) + x_s(:,t)*x_s(:,t)');
            end
    
            % Stato iniziale
            Sx1 = Sx1 + x_s(:,1);
            Sx1x1 = Sx1x1 + (P_s(:,:,1) + x_s(:,1)*x_s(:,1)');
    
            % Dinamica
            for t=2:Tn
                Sxxt   = Sxxt   + (P_s(:,:,t) + x_s(:,t)*x_s(:,t)');
                Sxxtm1 = Sxxtm1 + (P_sm_cross(:,:,t-1) + x_s(:,t)*x_s(:,t-1)');
                Sxtm1tm1 = Sxtm1tm1 + (P_s(:,:,t-1) + x_s(:,t-1)*x_s(:,t-1)');
            end
        end
    
        % --- M-step ---
        C = Syx / Sxx;
        d = (sum(cellfun(@(Y) size(Y,2), yCell)) > 0) * ...
            ((sum(cell2mat(cellfun(@(Y) sum(Y,2), yCell, 'uni', 0)),2) - C*sumX(yCell, xDim, C, d)) / Ttot);
        % Più semplice/robusto: ricalcola d dopo C:
        d = (1/Ttot) * (sum(cell2mat(cellfun(@(Y) sum(Y,2), yCell, 'uni', 0)),2) - C*sumXhatAcrossTrials(yCell, xDim, C, d)); %#ok<NASGU> % (puoi sostituire con accumulo esplicito nei loop)
    
        R = (Syy - C*Syx' - Syx*C' + C*Sxx*C') / Ttot;
    
        A = Sxxtm1 / Sxtm1tm1;
        Q = (Sxxt - A*Sxxtm1' - Sxxtm1*A' + A*Sxtm1tm1*A') / max(Ttot- N, 1);
    
        mu_1 = Sx1 / N;
        V_1  = Sx1x1/N - mu_1*mu_1';
    
        % simmetrizza/regularizza
        R = (R+R')/2; Q = (Q+Q')/2; V_1=(V_1+V_1')/2;
    
        % Convergenza
        if ll - ll_prev < tol*abs(ll_prev), break; end
        ll_prev = ll;
    end
    
    % Output
    kf.A=A; kf.C=C; kf.Q=Q; kf.R=R; kf.d=d; kf.mu_1=mu_1; kf.V_1=V_1;
end
