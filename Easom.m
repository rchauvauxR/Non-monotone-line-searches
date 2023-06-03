function value = Easom(option,x,D)
%   Implement the Easom function 
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

        value = -cos(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2);
        
    elseif option == "derivative"
%       Compute the gradient.
%       Parameters:
%           - x: The point being evaluated.
%           - D: The data of the function.

        value = zeros(1,2);
        value(1,1) = sin(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2) + cos(x(1))*cos(x(2))*2*(x(1)-pi)*exp(-(x(1)-pi)^2-(x(2)-pi)^2);
        value(1,2) = cos(x(1))*sin(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2) + cos(x(1))*cos(x(2))*2*(x(2)-pi)*exp(-(x(1)-pi)^2-(x(2)-pi)^2);
        
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