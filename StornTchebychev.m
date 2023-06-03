function value = StornTchebychev(option,x,D)
%   Implement the Storn's Tchebychev function 
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
        d = D{1};
        m = D{2};
        value = 0;
        u = 0;
        v = 0;
        for i = 1:n
            u = u + (1.2)^(n-i)*x(i);
            v = v + (-1.2)^(n-i)*x(i);
        end
        for j = 1:m
            ele = 0;
            for i = 1:n
                ele = ele + (2*j/m - 1)^(n-i)*x(i);
            end
            if ele > 1
                value = value + (ele - 1)^2;
            elseif ele < -1
                value = value + (ele + 1)^2;
            end
        end
        if u < d
            value = value + (u-d)^2;
        end
        if v < d
            value = value + (v-d)^2;
        end
        
    elseif option == "derivative"
%       Compute the gradient.
%       Parameters:
%           - x: The point being evaluated.
%           - D: The data of the function.

        [~,n] = size(x);
        d = D{1};
        m = D{2};
        value = zeros(1,n);
        u = 0;
        v = 0;
        for i = 1:n
            u = u + (1.2)^(n-i)*x(i);
            v = v + (-1.2)^(n-i)*x(i);
        end
        if u < d
            for i = 1:n
                value(1,i) = value(1,i) + 2*(u-d)*(1.2)^(n-i);
            end
        end
        if v < d
            for i = 1:n
                value(1,i) = value(1,i) + 2*(v-d)*(-1.2)^(n-i);
            end
        end
        
        for j = 1:m
            ele = 0;
            for i = 1:n
                ele = ele + (2*j/m - 1)^(n-i)*x(i);
            end
            if ele > 1
                for k = 1:n
                    value(1,k) = value(1,k) + 2*(ele-1)*(2*j/m-1)^(n-k);
                end
            elseif ele < -1
                for k = 1:n
                    value(1,k) = value(1,k) + 2*(ele+1)*(2*j/m-1)^(n-k);
                end
            end
        end
        
    elseif option == "data"
%       Provide the data of the function.
%       Parameters:
%           - x: Not applicable (no meaning).
%           - D: Not applicable (no meaning).

        d = 72.661;
        m = 60;
        value = {d,m};
        
    elseif option == "initial"
%       Provide a set of initial points.
%       Parameters:
%           - x: The number of desired initial points.
%           - D: The dimension of the initial points.
        
        % Grid 
        p = x/(2*D);
        a = -128;
        b = 128;
        
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