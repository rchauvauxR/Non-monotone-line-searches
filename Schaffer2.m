function value = Schaffer2(option,x,D)
%   Implement the Schaffer 2 function 
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

        value = (x(1)^2  + x(2)^2)^(0.25)*((sin(50*(x(1)^2+x(2)^2)^(0.1)))^2+1);
        
    elseif option == "derivative"
%       Compute the gradient.
%       Parameters:
%           - x: The point being evaluated.
%           - D: The data of the function.

        value = zeros(1,2);
        square = x(1)^2 + x(2)^2;
        sum1 = 0.5*(square)^(-0.75)*((sin(50*(square)^(0.1)))^2+1);
        sum2 = 20*square^(-0.65)*sin(50*square^(0.1))*cos(50*square^(0.1));
        ele = sum1 + sum2;
        value(1,1) = x(1)*ele;
        value(1,2) = x(2)*ele;
        
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
        a = -100;
        b = 100;
        
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