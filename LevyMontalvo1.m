function value = LevyMontalvo1(option,x,D)
%   Implement the Levy and Montalvo 1 function 
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
        y1 = (1+ 0.25*(x(1)+1));
        value = (pi/n)*10*sin(pi*y1)^2;
        
        for i = 1:n-1
            y2 = (1+ 0.25*(x(i+1)+1));
            value = value + pi/n*(y1 - 1)^2*(1+10*sin(pi*y2)^2);
            y1 = y2;
        end
        value = value + pi/n*(y2 - 1)^2;
        
    elseif option == "derivative"
%       Compute the gradient.
%       Parameters:
%           - x: The point being evaluated.
%           - D: The data of the function.

        [~,n] = size(x);
        value = zeros(1,n);
        y1 = (1+ 0.25*(x(1)+1)); 
        y2 = (1+ 0.25*(x(2)+1)); 
        dyi = 1/4;
        value(1,1) = pi/n*(20*pi*sin(pi*y1)*cos(pi*y1)*dyi + 2*(y1-1)*dyi*(1+10*sin(pi*y2)^2));
        for i = 2:n-1  
            y3 = (1+ 0.25*(x(i+1)+1));
            value(1,i) = pi/n*(2*(y2-1)*dyi*(1+10*sin(pi*y3)^2) + 20*pi*(y1-1)^2*sin(pi*y2)*cos(pi*y2)*dyi);
            y1 = y2;
            y2 = y3;
        end
        value(1,n) = pi/n*(2*(y2-1)*dyi + 20*pi*(y1-1)^2*sin(pi*y2)*cos(pi*y2)*dyi);
        
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


