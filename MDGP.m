function value = MDGP(option,x,D)
%   Implement the Molecular Distance Geometry Problem function 
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

        [M,~] = size(D{1});
        [~,N] = size(x);
        N = N/3;
        value = 0;
        vec = zeros(1,3);
        for i=1:M
            v1 = D{1}(i,1);
            v2 = D{1}(i,2);         
            vec(1) = x(3*v1 -2) - x(3*v2 -2);
            vec(2) = x(3*v1 -1) - x(3*v2 -1);
            vec(3) = x(3*v1) - x(3*v2);
            value = value + ((norm(vec))^2 - D{1}(i,3)^2)^2;
        end

    elseif option == "derivative"
%       Compute the gradient.
%       Parameters:
%           - x: The point being evaluated.
%           - D: The data of the function.

        data = D{1};
        before = D{2};
        after = D{3};
        distance = D{4};
        [~,N] = size(x);
        N = N/3;
        vec = zeros(3,1);
        gradient=zeros(1,3*N);
        for i = 1:N
            Comp1 = 0;
            Comp2 = 0;
            Comp3 = 0;
            for j = after{i}  
                vec(1) = x(3*i -2) - x(3*j -2);
                vec(2) = x(3*i -1) - x(3*j -1);
                vec(3) = x(3*i) - x(3*j);
                Comp1 = Comp1 + 4*(x(3*i -2) - x(3*j -2))*((norm(vec))^2 - distance(i,j)^2);
                Comp2 = Comp2 + 4*(x(3*i -1) - x(3*j -1))*((norm(vec))^2 - distance(i,j)^2);
                Comp3 = Comp3 + 4*(x(3*i) - x(3*j))*((norm(vec))^2 - distance(i,j)^2);             
            end
            for j = before{i}
                vec(1) = x(3*j -2) - x(3*i -2);
                vec(2) = x(3*j -1) - x(3*i -1);
                vec(3) = x(3*j) - x(3*i);
                Comp1 = Comp1 - 4*(x(3*j -2) - x(3*i -2))*((norm(vec))^2 - distance(j,i)^2);
                Comp2 = Comp2 - 4*(x(3*j -1) - x(3*i -1))*((norm(vec))^2 - distance(j,i)^2);
                Comp3 = Comp3 - 4*(x(3*j) - x(3*i))*((norm(vec))^2 - distance(j,i)^2);             
            end
            gradient(3*i -2) = Comp1;
            gradient(3*i -1) = Comp2;
            gradient(3*i) = Comp3;
        end     
        value = gradient;
        
    elseif option == "data"
%       Provide the data of the function.
%       Parameters:
%           - x: Name of a file containg the protein distances.
%           - D: Not applicable (no meaning).
        
        Data = fopen(x);
        line = fgetl(Data);
        liste = strsplit(line);
        N = str2double(liste(1));
        line = fgetl(Data);
        liste = strsplit(line);
        M = str2double(liste(1));
        line = fgetl(Data);
        liste = strsplit(line);
        min = str2double(liste(1));
        line = fgetl(Data);
        liste = strsplit(line);
        max = str2double(liste(1));
        Matrix = zeros(M,3);

        i = 1;
        while ~feof(Data)
            line = fgetl(Data);
            liste = strsplit(line);
            Matrix(i,1)=str2double(liste(1));
            Matrix(i,2)=str2double(liste(2));
            Matrix(i,3)=str2double(liste(3));
            i=i+1;
        end
        
        before = cell(N,1);
        after = cell(N,1);
        distance = zeros(N);
        for e = 1:M
            v1 = int64(Matrix(e,1));
            v2 = int64(Matrix(e,2));
            after{v1} = [after{v1}, v2];
            before{v2} = [before{v2}, v1];
            distance(v1,v2) = Matrix(e,3);
            distance(v2,v1) = Matrix(e,3);
        end
        value = {Matrix,before,after,distance,min,max};
        fclose(Data);
        
    elseif option == "initial"
%       Provide a set of initial points.
%       Parameters:
%           - x: The number of desired initial points.
%           - D: The data of the function.
      
        rng(1,'twister');
        [N,~] = size(D{2});
        a = D{5};
        b = D{6};
        for i = 1:x
            x0(i,:) = a + (b-a)*rand(1,3*N);
        end
        value = x0;       

    else
        printf("error");
    end
end