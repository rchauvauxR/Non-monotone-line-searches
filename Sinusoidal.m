function value = Sinusoidal(option,x,D)
%   Implement the Sinusoidal function 
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
        A = D(1);
        B = D(2);
        z = D(3);
        prod1 = 1;
        prod2 = 1;
        for i = 1:n
            prod1 = prod1*sin(x(i)-z);
            prod2 = prod2*sin(B*(x(i)-z));
        end
        value = -(A*prod1 + prod2);
        
    elseif option == "derivative"
%       Compute the gradient.
%       Parameters:
%           - x: The point being evaluated.
%           - D: The data of the function.

        [~,n] = size(x);
        A = D(1);
        B = D(2);
        z = D(3);
        value = ones(1,n);
        for i = 1:n
            prod1 = 1;
            prod2 = 1;
            for k = 1:n
                if k == i
                    prod1 = prod1*cos(x(i) - z);
                    prod2 = prod2*B*cos(B*(x(i)-z));
                else
                    prod1 = prod1*sin(x(k) - z);
                    prod2 = prod2*sin(B*(x(k)-z));
                end
            end
            value(1,i) = -(A*prod1 + prod2);
        end
        
    elseif option == "data"
%       Provide the data of the function.
%       Parameters:
%           - x: Not applicable (no meaning).
%           - D: Not applicable (no meaning).

         value = [2.5,5,30];
         
    elseif option == "initial"
%       Provide a set of initial points.
%       Parameters:
%           - x: The number of desired initial points.
%           - D: The dimension of the initial points.

        % Grid
        p = x/(2*D);
        a = 0;
        b = 180;
        
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
