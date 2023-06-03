function value = PriceTransistorModelling(option,x,D)
%   Implement the Price's Transistor Modelling function 
%   This function provides four possibilities based on the option parameter:
%       - "function"   : Compute the function value
%       - "derivative" : Compute the gradient
%       - "data"       : Provide the data of the function
%       - "initial"    : Provide a set of initial points
%   The interpretation of the parameters x and D depends on the chosen option.
%   The function returns value, which represents the desired result based on the chosen option.

    if option == "function"
%       Compute the function value.
%       Parameters:
%           - x: The point being evaluated.
%           - D: The data of the function.

        g = D;
        gamma = x(1)*x(3) - x(2)*x(4);
        value = gamma^2;
        for k = 1:4
            alphak = (1-x(1)*x(2))*x(3)*(exp(x(5)*(g(1,k) - g(3,k)*x(7)*10^(-3) - g(5,k)*x(8)*10^(-3)))-1) - g(5,k) + g(4,k)*x(2);
            betak = (1-x(1)*x(2))*x(4)*(exp(x(6)*(g(1,k) - g(2,k) - g(3,k)*x(7)*10^(-3) + g(4,k)*x(9)*10^(-3)))-1) - g(5,k)*x(1) + g(4,k);
            value = value + alphak^2 + betak^2;
        end
        
    elseif option == "derivative"
%       Compute the gradient.
%       Parameters:
%           - x: The point being evaluated.
%           - D: The data of the function.

        g = D;
        value = zeros(1,9);
        gamma = x(1)*x(3) - x(2)*x(4);
        value(1,1) = 2*gamma*x(3);
        value(1,2) = -2*gamma*x(4);
        value(1,3) = 2*gamma*x(1);
        value(1,4) = -2*gamma*x(2);
        for k = 1:4
            alphak = (1-x(1)*x(2))*x(3)*(exp(x(5)*(g(1,k) - g(3,k)*x(7)*10^(-3) - g(5,k)*x(8)*10^(-3)))-1) - g(5,k) + g(4,k)*x(2);
            betak = (1-x(1)*x(2))*x(4)*(exp(x(6)*(g(1,k) - g(2,k) - g(3,k)*x(7)*10^(-3) + g(4,k)*x(9)*10^(-3)))-1) - g(5,k)*x(1) + g(4,k);
            
            value(1,1) = value(1,1) - 2*alphak*x(2)*x(3)*(exp(x(5)*(g(1,k) - g(3,k)*x(7)*10^(-3) - g(5,k)*x(8)*10^(-3)))-1) - 2*betak*(x(2)*x(4)*(exp(x(6)*(g(1,k) - g(2,k) - g(3,k)*x(7)*10^(-3) + g(4,k)*x(9)*10^(-3)))-1) + g(5,k));
            value(1,2) = value(1,2) - 2*alphak*(x(1)*x(3)*(exp(x(5)*(g(1,k) - g(3,k)*x(7)*10^(-3) - g(5,k)*x(8)*10^(-3)))-1) - g(4,k)) - 2*betak*x(1)*x(4)*(exp(x(6)*(g(1,k) - g(2,k) - g(3,k)*x(7)*10^(-3) + g(4,k)*x(9)*10^(-3)))-1);
            value(1,3) = value(1,3) + 2*alphak*(1-x(1)*x(2))*(exp(x(5)*(g(1,k) - g(3,k)*x(7)*10^(-3) - g(5,k)*x(8)*10^(-3)))-1);
            value(1,4) = value(1,4) + 2*betak*(1-x(1)*x(2))*(exp(x(6)*(g(1,k) - g(2,k) - g(3,k)*x(7)*10^(-3) + g(4,k)*x(9)*10^(-3)))-1);
            value(1,5) = value(1,5) + 2*alphak*(1-x(1)*x(2))*x(3)*(g(1,k) - g(3,k)*x(7)*10^(-3) - g(5,k)*x(8)*10^(-3))*exp(x(5)*(g(1,k) - g(3,k)*x(7)*10^(-3) - g(5,k)*x(8)*10^(-3)));
            value(1,6) = value(1,6) + 2*betak*(1-x(1)*x(2))*x(4)*(g(1,k) - g(2,k) - g(3,k)*x(7)*10^(-3) + g(4,k)*x(9)*10^(-3))*exp(x(6)*(g(1,k) - g(2,k) - g(3,k)*x(7)*10^(-3) + g(4,k)*x(9)*10^(-3)));
            value(1,7) = value(1,7) - 2*alphak*(1-x(1)*x(2))*x(3)*g(3,k)*x(5)*10^(-3)*exp(x(5)*(g(1,k) - g(3,k)*x(7)*10^(-3) - g(5,k)*x(8)*10^(-3))) -2*betak*(1-x(1)*x(2))*x(4)*g(3,k)*x(6)*10^(-3)*exp(x(6)*(g(1,k) - g(2,k) - g(3,k)*x(7)*10^(-3) + g(4,k)*x(9)*10^(-3)));
            value(1,8) = value(1,8) - 2*alphak*(1-x(1)*x(2))*x(3)*g(5,k)*x(5)*10^(-3)*exp(x(5)*(g(1,k) - g(3,k)*x(7)*10^(-3) - g(5,k)*x(8)*10^(-3)));
            value(1,9) = value(1,9) + 2*betak*(1-x(1)*x(2))*x(4)*g(4,k)*x(6)*10^(-3)*exp(x(6)*(g(1,k) - g(2,k) - g(3,k)*x(7)*10^(-3) + g(4,k)*x(9)*10^(-3)));     
        end
        
    elseif option == "data"
%       Provide the data of the function.
%       Parameters:
%           - x: Not applicable (no meaning).
%           - D: Not applicable (no meaning).

        value = [[0.485 0.752 0.869 0.982],
                [0.369 1.254 0.703 1.455],
                [5.2095 10.0677 22.9274 20.2153],
                [23.3037 101.779 111.461 191.267],
                [28.5132 111.8467 134.3884 211.4823]];
            
    elseif option == "initial"
%       Provide a set of initial points.
%       Parameters:
%           - x: The number of desired initial points.
%           - D: The dimension of the initial points.

        % Grid
        p = x/(2*D);
        a = -10;
        b = 10;
        
        value = zeros(x,D);
        index = 1;
        x_mid = zeros(1,D) + (b-a)/2;
        for j = 1:p
            for k = 1:D
                unit = zeros(1,D);
                unit(1,k)= 1;
                value(index,:) = x_mid + unit*j*(b-a)/(2*p);
                index = index+1;
                value(index,:) = x_mid - unit*j*(b-a)/(2*p);
                index = index+1;
            end
        end
        
        % Random
        %{
        value = zeros(x,D);
        for i = 1:x
            value(i,:) = a + (b-a)*rand(1,D);
        end
        %}
        
    else
        printf("error");
    end
end