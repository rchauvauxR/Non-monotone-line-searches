function value = LevyMontalvo2(option,x,D)
%   Implement the Levy and Montalvo 2 function 
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
        value = 0.1*(sin(3*pi*x(1))^2 + (x(n) - 1)^2*(1 + sin(2*pi*x(n))^2));
        
        for i = 1:n-1
            value = value + 0.1*(x(i) - 1)^2*(1+sin(3*pi*x(i+1))^2);
        end
        
    elseif option == "derivative"
%       Compute the gradient.
%       Parameters:
%           - x: The point being evaluated.
%           - D: The data of the function.

        [~,n] = size(x);
        value = zeros(1,n);
        value(1,1) = 0.1*(6*pi*sin(3*pi*x(1))*cos(3*pi*x(1)) + 2*(x(1)-1)*(1+sin(3*pi*x(2))^2));
        for i = 2:n-1  
            value(1,i) = 0.1*(2*(x(i)-1)*(1+sin(3*pi*x(i+1))^2) + 6*pi*(x(i-1)-1)^2*sin(3*pi*x(i))*cos(3*pi*x(i)));
        end
        value(1,n) = 0.1*(2*(x(n)-1)*(1+sin(2*pi*x(n))^2) + 4*pi*(x(n)-1)^2*sin(2*pi*x(n))*cos(2*pi*x(n))+ 6*pi*(x(n-1)-1)^2*sin(3*pi*x(n))*cos(3*pi*x(n)));
        
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
        a = -5;
        b = 5;
        
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


