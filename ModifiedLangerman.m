function value = ModifiedLangerman(option,x,D)
%   Implement the Modified Langerman function 
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
        c = D{1};
        a = D{2};
        value = 0;
        for j = 1:5
            dj = 0;
            for i = 1:n
                dj = dj + (x(i) - a(j,i))^2;
            end
            value = value - c(j)*cos(dj/pi)*exp(-pi*dj);
        end
        
    elseif option == "derivative"
%       Compute the gradient.
%       Parameters:
%           - x: The point being evaluated.
%           - D: The data of the function.

        [l,n] = size(x);
        c = D{1};
        a = D{2};
        value = zeros(1,n);
        d = zeros(5,1);
        for j = 1:5
            for i = 1:n
                d(j,1) = d(j,1) + (x(i) - a(j,i))^2;
            end
        end
        for i = 1:n
            for j = 1:5
                dd=2*(x(i) - a(j,i));
                value(1,i) = value(1,i) + c(j)*(sin(d(j,1)/pi)/pi*dd*exp(-pi*d(j,1)) + cos(d(j,1)/pi)*pi*dd*exp(-pi*d(j,1)));
            end     
        end
        
    elseif option == "data"
%       Provide the data of the function.
%       Parameters:
%           - x: Not applicable (no meaning).
%           - D: Not applicable (no meaning).

        a = [[9.681 0.667 4.783 9.095 3.517 9.325 6.544 0.211 5.122 2.020],
             [9.400 2.041 3.788 7.931 2.882 2.672 3.568 1.284 7.033 7.374],
             [8.025 9.152 5.114 7.621 4.564 4.711 2.996 6.126 0.734 4.982],
             [2.196 0.415 5.649 6.979 9.510 9.166 6.304 6.054 9.377 1.426],
             [8.074 8.777 3.467 1.867 6.708 6.349 4.534 0.276 7.633 1.567]];
        c = [0.806 0.517 0.100 0.908 0.965];
        value = {c,a};
        
    elseif option == "initial"
%       Provide a set of initial points.
%       Parameters:
%           - x: The number of desired initial points.
%           - D: The dimension of the initial points.

        % Grid
        p = x/(2*D);
        a = 0;
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