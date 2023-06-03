function value = EpistaticMichalewicz(option,x,D)
%   Implement the Epistatic Michalewicz function 
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

        m = D(1);
        theta = D(2);
        [l,n] = size(x);
        value = 0;
        for i = 1:n
            if i == n
                y = x(n);
            elseif mod(i,2) == 0
                y = x(i)*sin(theta) + x(i+1)*cos(theta);
            else
                y = x(i)*cos(theta) - x(i+1)*sin(theta);
            end
            value = value - sin(y)*(sin(i*y^2/pi))^(2*m);
        end
        
    elseif option == "derivative"
%       Compute the gradient.
%       Parameters:
%           - x: The point being evaluated.
%           - D: The data of the function.

        m = D(1);
        theta = D(2);
        [l,n] = size(x);
        value = zeros(1,n);
        for i = 1:n
            if i == n
                y = x(n);
                dy = 1;
                if mod(n,2)==0 
                    dy1 = -sin(theta);
                    y1 = x(i-1)*cos(theta) - x(i)*sin(theta);
                else
                    dy1 = cos(theta);
                    y1 = x(i-1)*sin(theta) + x(i)*cos(theta);
                end
            elseif mod(i,2) == 0
                y = x(i)*sin(theta) + x(i+1)*cos(theta);
                dy = sin(theta);
                dy1 = -sin(theta);
                y1 = x(i-1)*cos(theta) - x(i)*sin(theta);
            else
                y = x(i)*cos(theta) - x(i+1)*sin(theta);
                dy = cos(theta);
                if i == 1
                    dy1 = 0;
                    y1 = 0;
                else
                    dy1 = cos(theta);
                    y1 = x(i-1)*sin(theta) + x(i)*cos(theta);
                end
            end
            part1 = dy*cos(y)*(sin(i*y^2/pi))^(2*m);
            part2 = 4*m*i*y/pi*dy*sin(y)*cos(i*y^2/pi)*(sin(i*y^2/pi))^(2*m-1);
            part3 = dy1*cos(y1)*(sin((i-1)*y1^2/pi))^(2*m);
            part4 = 4*m*(i-1)*y1/pi*dy1*sin(y1)*cos((i-1)*y1^2/pi)*(sin((i-1)*y1^2/pi))^(2*m-1);
            value(1,i) = -part1 - part2 - part3 - part4;
        end
        
    elseif option == "data"
%       Provide the data of the function.
%       Parameters:
%           - x: Not applicable (no meaning).
%           - D: Not applicable (no meaning).

        value = [10,pi/6];
        
    elseif option == "initial"
%       Provide a set of initial points.
%       Parameters:
%           - x: The number of desired initial points.
%           - D: The dimension of the initial points.

        % Grid
        p = x/(2*D);
        a = 0;
        b = pi;
        
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


