function [Historic,ind_hist] = NM2_alg(MDGP,Problem,D,x0, alpha0, beta, rho,f_min,Param)
% Implement the non-monotone Algorithm NM2.
% 
% Parameters:
%   - MDGP      : Set to 0 if it is not used on an MDGP problem.
%   - Problem   : The problem being solved.
%   - D         : The data of the problem.
%   - x0        : An array of doubles representing the initial point.
%   - alpha0    : A positive double.
%   - beta      : A double in the range (0,1).
%   - rho       : A double in the range (0,1).
%   - f_min     : The optimal function value of the problem.
%   - Param     : Not applicable (no meaning).
%
% Outputs:
%   - Historic  : An array of doubles representing the found function values.
%   - ind_hist  : The number of function evaluations used.
    
    % Step 0
    [~,np] = size(x0);
    kp_max = 100*(np+1);
    if MDGP == 0
        Historic = zeros(1,1100) + Inf;
    else
        Historic = zeros(1,kp_max) + Inf;
    end
    ind_hist = 1;
    
    dim = length(x0);
    I = eye(dim);
    Hk = I;
    
    k = 0; 
    xk = x0;
    xk1 = xk;
    
    alphak = alpha0;
    gk = Problem("derivative",x0,D);
    gk1 = gk;
    
    C_k = 0;
    Q_k = 1;
    
    tol = 10^(-300);
    err =1;
    
    k_check = 0;
    f_check = 0;

    while(ind_hist < kp_max && err > tol )
        
        % Step 1
        if k == 0
            f_k = Problem("function",xk1,D);
            Historic(1,ind_hist) = f_k;
            ind_hist = ind_hist +1;
        else
            f_k = f_kl;
        end
        
        if MDGP ~= 0 
            if f_k <= f_min
                break;
            end
            if ind_hist - k_check > np/10
                if abs(f_check - f_k) < 10^(-1)
                    break;
                else
                    k_check = ind_hist;
                    f_check = f_k;
                end
            end
        end
        % BFGS update
        if k ~= 0
            sk = xk1 - xk;
            yk = gk1 - gk;
            prod = sk*yk';
            if prod > 0
                Hk = (I - sk'*yk/prod)*Hk*(I - yk'*sk/prod) + sk'*sk/prod;
            end
            Q_k1 = (0.85/k)*Q_k + 1;
            C_k = ((0.85/k)*Q_k*C_k + f_k)/Q_k1;
            Q_k = Q_k1;
            
        end
        gk = gk1;
        xk = xk1;
        dk = -Hk*gk';
        
        % Step 2.1
        Notfind = 1;
        l = 0;
        v_k = max(C_k - f_k,0);
        
        % Step 2.2
        while(Notfind)
            if ind_hist == kp_max 
                break
            end
            f_kl = Problem("function",xk + alphak*beta^l*dk',D);
            Historic(1,ind_hist) = f_k;
            ind_hist = ind_hist +1;
            if f_kl <= f_k + rho*alphak*beta^l*gk*dk + v_k
                break;
            else
                l = l+1;
            end
        end
        
        % Step 3
        xk1 = xk + alphak*beta^l*dk';
        alphak = alphak*beta^(l-1);
        gk1 = Problem("derivative",xk1,D);
        k = k+1;  
        err = norm(gk);
    end
end