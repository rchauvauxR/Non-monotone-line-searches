function value = Neumaier2(option,x,D)
%   Implement the Neumaier 2 function 
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

        [~,n] = size(x);
        b = D{1};
        value = 0;
        for k = 1:n
            sum = 0;
            for i = 1:n
                sum = sum + x(i)^k;
            end
            value = value + (b(k) - sum)^2;
        end
        
    elseif option == "derivative"
%       Compute the gradient.
%       Parameters:
%           - x: The point being evaluated.
%           - D: The data of the function.

        [~,n] = size(x);
        value = zeros(1,n);
        b = D{1};
        for k = 1:n
            sum = 0;
            for i = 1:n
                sum = sum + x(i)^k;
            end
            ele = -2*k*(b(k) - sum);
            for i = 1:n
                value(1,i) = value(1,i) + ele*x(i)^(k-1);
            end
        end
        
    elseif option == "data"
%       Provide the data of the function.
%       Parameters:
%           - x: Not applicable (no meaning).
%           - D: Not applicable (no meaning).
        b = [8,18,44,114];
        value = {b};
        
    elseif option == "initial"
%       Provide a set of initial points.
%       Parameters:
%           - x: The number of desired initial points.
%           - D: The dimension of the initial points.

        % Grid
        
        p = x/(2*D);
        a = 0;
        b = D;
        
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