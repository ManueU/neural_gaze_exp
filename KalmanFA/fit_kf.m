function params = fit_kf(X, z, ridge_lambda)
        [T, d] = size(X);
        
        Xprev = X(1:end-1, :);
        Xnext = X(2:end,   :);
        rowsOK = all(isfinite([Xprev Xnext]), 2);
        Xprev = Xprev(rowsOK, :);
        Xnext = Xnext(rowsOK, :);
        mu = mean(Xprev,1);
        Xp = Xprev - mu;
        Xn = Xnext - mu;

        G = Xp'*Xp + ridge_lambda*eye(size(X,2));
        A  = (G \ (Xp' * Xn))';
        b  = (eye(size(A))-A)*mu';                                   
        
        Ew = Xnext - (Xprev*A' + b');
        Q  = cov(Ew,1) + 1e-9*eye(size(Ew,2));                          
        
        X1 = [X, ones(T,1)];    
        if ridge_lambda > 0
            I = eye(d+1); I(end,end) = 0;          
            W = (X1' * X1 + ridge_lambda * I) \ (X1' * z); 
        else
            W = X1 \ z;                             
        end
        C = W(1:d, :)';                             
        dvec = W(end, :)';                         
        Ez = z - (X * C' + dvec');                
        R  = cov(Ez, 1);                            
        
        x0 = mean(X(1:10, :), 1)';                  
        P0 = eye(d) * 1e-2;                         
        
        params = struct('A',A,'Q',Q,'C',C,'d',dvec,'R',R,'x0',x0,'P0',P0);
end
