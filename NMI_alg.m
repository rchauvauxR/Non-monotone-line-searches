function [Historic,ind_hist] = NMI_alg(MDGP,Problem,D,x0, alpha0, beta, rho,f_min,Param)
% Implement the non-monotone informed Algorithm NMI.
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
%   - Param     : [r,M,delta]
%
% Outputs:
%   - Historic  : An array of doubles representing the found function values.
%   - ind_hist  : The number of function evaluations used.  

    %Step 0 
    
    r = Param(1);
    M = Param(2);
    d = Param(3);
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
    f_k = Problem("function",xk,D);            
    f_0 = f_k;
    Historic(1,ind_hist) = f_k;
    ind_hist = ind_hist +1;

    alphak = alpha0;
    gk = Problem("derivative",x0,D);
    
    R = r;   
    phi = 1.01;
    sigma = M/d^2;
    err = 1;
    tol = 10^(-300);
    v_k = 0;
    
    while(ind_hist < kp_max && err > tol)
        
        
        % Step 1
        dk = -Hk*gk';
        
        % Step 2.1
        Notfind = 1;
        l = 0;
        
        % Step 2.2
        while(Notfind)
            if ind_hist == kp_max 
                break
            end
            xk2 = xk + alphak*beta^l*dk';         
            f_kl = Problem("function",xk2,D);
            Historic(1,ind_hist) = f_kl;
            ind_hist = ind_hist +1;
            if f_kl <= f_k + rho*alphak*beta^l*gk*dk + v_k
                break;
            else
                l = l+1;
            end
        end
        
        
        % Step 3.1
        f_k = f_kl;
        xk1 = xk2;
        gk1 = Problem("derivative",xk1,D);
        
        % Step 3.2
        c = 0;
        %if f_k-f_min > d^2*(f_0 - f_min) && norm(gk1) <d*max(1,(f_k-f_min))
        if f_k-f_min > d^2*(f_0 - f_min) && norm(gk1) <d*min(f_0 - f_min,f_k-f_min)
            c = 1;
        end
        
        % Step 3.3-3.4
        if c == 1
            alphak = R/norm(gk1);
            v_k = sigma*min(f_0 - f_min,f_k-f_min)*(1/(1+k))^(phi);
            Hk = I;
        else
            alphak = alphak*beta^(l-1);
            v_k = 0;
            sk = xk1 - xk;
            yk = gk1 - gk;
            prod = sk*yk';
            if prod > 0
                Hk = (I - sk'*yk/prod)*Hk*(I - yk'*sk/prod) + sk'*sk/prod;
            end         
        end
        
        % Step 4
        if MDGP ~= 0
            if f_k <= f_min
                break;
            end
        end
        gk = gk1;
        xk = xk1;
        k = k+1;  
        err = norm(gk);
    end
end