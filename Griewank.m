function value = Griewank(option,x,D)
%   Implement the Griewank function 
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

        [l,n] = size(x);
        value = 1;
        sum = 0;
        product = 1;
        for i = 1:n
            sum = sum + x(i)^2;
            product = product*cos(x(i)/sqrt(i));
        end
        value = value + sum/4000 - product;

    elseif option == "derivative"
%       Compute the gradient.
%       Parameters:
%           - x: The point being evaluated.
%           - D: The data of the function.

        [l,n] = size(x);
        value = zeros(1,n);
        for i = 1:n
            product = 1;
            for j = 1:n
                if j ~= i
                    product = product*cos(x(j)/sqrt(j));
                end
            end
            value(1,i) = x(i)/2000 + sin(x(i)/sqrt(i))/sqrt(i)*product;
        end
        
    elseif option == "data"
%       Provide the data of the function.
%       Parameters:
%           - x: Not applicable (no meaning).
%           - D: Not applicable (no meaning).

        value = 0;
        
    elseif option == "initial"
%       Provide a set of initial points.
%       Parameters:
%           - x: The number of desired initial points.
%           - D: The dimension of the initial points.

        % Grid
        p = x/(2*D);
        a = -600;
        b = 600;
        
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
        
        % Initial Grid
        %{
        value = zeros(60,2);
        index = 1;
        for i = 1:4
            for j = 1:15
                value(index,1) = -600 + 1200*(i-1)/3;
                value(index,2) = -600 + 1200*(j-1)/14;
                index = index +1;
            end
        end
        %}
        
    else
        printf("error");
    end
end