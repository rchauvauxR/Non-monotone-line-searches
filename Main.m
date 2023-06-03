% Author: Romain Chauvaux
% Date: 27/05/2023
% 
% This main file allows the user to conduct numerical experiments.
% 
% Instructions:
% 1. Choose the methods you want to compare.
% 2. Select the problems you want to use. Note that the MDGP problem cannot be
%    mixed with other problems.
% 3. Call the main function with your desired choices.



% Choice of the methods and their Param

Methods = {@M1_alg, @NM1_alg, @NM2_alg, @NM3_alg, @NMI_alg};
Names = ["M1","NM1","NM2","NM3","NMI"];
Param_M1 = [];
Param_NM1 = [];
Param_NM2 = [];
Param_NM3 = [0,2.1135]; % [Sigma, Theta] remark : Set Sigma = 0 for Sigma = |f(x_0)|
Param_NMI = [1,100,0.001]; % [R, M, delta]
%Param_NMI = [10,10^9,0.001];
Param = {Param_M1,Param_NM1,Param_NM2,Param_NM3,Param_NMI};

% Choice of the desired test problems 
B1 = {@Bohachevsky1,2,0};
B2 = {@Bohachevsky2,2,0};
CM = {@CosineMixture,4,-0.4};
EP = {@Easom,2,-1};
EM = {@EpistaticMichalewicz,10,-9.660152};
EXP = {@Exponential,10,-1};
GW = {@Griewank,2,0};
LM1 = {@LevyMontalvo1,3,0};
LM2 = {@LevyMontalvo2,10,0};
ML = {@ModifiedLangerman,10,-0.965};
NF2 = {@Neumaier2,4,0};
NF3 = {@Neumaier3,10,-210.000000001};
PTM = {@PriceTransistorModelling,9,0};
RG = {@Rastrigin,10,0};
SF1 = {@Schaffer1,2,0};
SF2 = {@Schaffer2,2,0}; 
FX = {@ShekelsFoxholes,10,-10.208792792153845};
SBT = {@Shubert,2,-186.7309088310238};
SIN = {@Sinusoidal,10,-3.5};
ST = {@StornTchebychev,9,0};


R = 360;
Problems = {B1,B2,CM,EP,EM,EXP,GW,LM1,LM2,ML,NF2,NF3,PTM,RG,SF1,SF2,FX,SBT,SIN,ST};
Names_prob = ["Bohachevsky1","Bohachevsky2","CosineMixture","Easom","EpistaticMichalewicz","Exponential","Griewank","LevyMontalvo1","LevyMontalvo2","ModifiedLangerman","Neumaier2","Neumaier3","PriceTransistorModelling","Rastrigin","Schaffer1","Schaffer2","ShekelsFoxholes","Shubert","Sinusoidal","StornTchebychev"];
[~,N] = size(Problems);

Data = cell(1,N);
x0 = cell(N,R);
for i = 1:N
    D = Problems{i}{1}("data",1,1);
    [l1,l2] = size(D);
    Data(i) = mat2cell(D,l1,l2);
    x0(i,:) = mat2cell(Problems{i}{1}("initial",R,Problems{i}{2}),ones(R,1),Problems{i}{2});
end

%main(Methods,Names,Param,Problems,Names_prob,Data,x0,"Box");
%main(Methods,Names,Param,Problems,Names_prob,Data,x0,"SA");
main(Methods,Names,Param,Problems,Names_prob,Data,x0,"Profiles");

% Choice of the MDGP 
R = 2;
Problems = {@MDGP};
Names_prob = ["MDGP"];
protein = "PythonCode/1ao2_5.txt";
Data = {Problems{1}("data",protein,0)};
[M,~] = size(Data{1}{1});
[N,~] = size(Data{1}{2});
x0 = cell(1,R);
x0(1,:) = mat2cell(Problems{1}("initial",R,Data{1}),ones(R,1),3*N);

%main(Methods,Names,Param,Problems,Names_prob,Data,x0,"Box_MDGP");
%main(Methods,Names,Param,Problems,Names_prob,Data,x0,"SA_MDGP");
%main(Methods,Names,Param,Problems,Names_prob,Data,x0,"Profiles_MDGP");



function main(Methods,Names,Param,Problems,Names_prob,Data,x0,choice)
%   The main fucntion
%   This function provides six possibilities based on the choice parameter:
%       - "Box"           : Compare each methods with a boxplot per problems.
%       - "Box_MDGP"      : Compare each methods with a boxplot on the MDGP.
%       - "SA"            : Compare each methods and Simulated Annealing with a boxplot per problems.
%       - "SA_MDGP"       : Compare each methods and Simulated Annealing with a boxplot on the MDGP.
%       - "Profiles"      : Compare each methods with a data and performance profiles.
%       - "Profiles_MDGP" : Compare each methods with a data and performance profiles on the MDGP.

    [N,R] = size(x0);
    alpha0 = 1;
    beta = 0.5;
    rho = 0.5;
    gate = 10^(-5);
    N2 = zeros(1,N*R);
    for i = 1:N
        [~,dim] = size(x0{i,1});
        N2(1,(i-1)*R+1:i*R) = ones(1,R)*dim + 1;
    end
    switch choice 
        case "Box"
            H = Box(Methods,Param,Problems,Data,x0,alpha0,beta,rho);
            Box_plot(H,Names,Names_prob);
        case "Box_MDGP"
            H = Box_MDGP(Methods,Param,Problems,Data,x0,alpha0,beta,rho);
            Box_plot(H,Names,Names_prob);
        case "SA"
            H = SA(Methods,Param,Problems,Data,x0,alpha0,beta,rho);
            Names = [Names,"SA"];
            Box_plot(H,Names,Names_prob);
        case "SA_MDGP"
            H = SA_MDGP(Methods,Param,Problems,Data,x0,alpha0,beta,rho);
            Names = [Names,"SA"];
            Box_plot(H,Names,Names_prob);
        case "Profiles"
            H = Profiles(Methods,Param,Problems,Data,x0,alpha0,beta,rho);
            Data_profile(H,N2,gate,Names); 
            perf_profile(H,gate,1,Names);  
        case "Profiles_MDGP"
            H = Profiles_MDGP(Methods,Param,Problems,Data,x0,alpha0,beta,rho);
            Data_profile(H,N2,gate,Names); 
            perf_profile(H,gate,1,Names); 
        otherwise
            disp("Parameter not recognised");
    end
            
    
end

function H = Box(Methods,Param,Problems,Data,x0,alpha0,beta,rho)
% The function Box
%
% Parameters:
%   - Methods  : The desired methods to compare.
%   - Param    : The specific parameters for each method.
%   - Problems : The considered problems.
%   - Data     : The specific data for each problem.
%   - x0       : The initial points for each test problem.
%   - alpha0   : A positive double.
%   - beta     : A double in the range (0,1).
%   - rho      : A double in the range (0,1).
%
% Outputs:
%   - H        : An array of doubles representing the minimal of the found function values for each method on each test problem.

    MDGP = 0;
    [N,R] = size(x0);
    [~,K] = size(Methods);
    H = zeros(N*R,K);
    for i = 1:N
        for j = 1:R
            index = (i-1)*R+j;
            for k = 1:K
                [List,~] = Methods{k}(MDGP,Problems{i}{1},Data{i},x0{i,j},alpha0,beta,rho,Problems{i}{3},Param{k});
                H(index,k) = min(List);
            end
        end
    end
end

function H = Box_MDGP(Methods,Param,Problems,Data,x0,alpha0,beta,rho)
% The function Box_MDGP
%
% Parameters:
%   - Methods  : The desired methods to compare.
%   - Param    : The specific parameters for each method.
%   - Problems : The MDGP problems.
%   - Data     : The specific data for MDGP.
%   - x0       : The initial points for each test problem.
%   - alpha0   : A positive double.
%   - beta     : A double in the range (0,1).
%   - rho      : A double in the range (0,1).
%
% Outputs:
%   - H        : An array of doubles representing the minimal of the found function values for each method on the MDGP.

    MDGP = 1;
    [M,~] = size(Data{1}{1});
    [N,R] = size(x0);
    [~,K] = size(Methods);
    
    epsilon = 10^(-5);
    alpha0 = 1/M;
    
    fmin = M*epsilon;
    H = zeros(N*R,K);
    for i = 1:N
        for j = 1:R
            index = (i-1)*R+j;
            for k = 1:K
                [List,~] = Methods{k}(MDGP,Problems{i},Data{i},x0{i,j},alpha0,beta,rho,fmin,Param{k});
                H(index,k) = min(List);
            end
        end
    end
end

function H = SA(Methods,Param,Problems,Data,x0,alpha0,beta,rho)
% The function SA
%
% Parameters:
%   - Methods  : The desired methods to compare.
%   - Param    : The specific parameters for each method.
%   - Problems : The considered problems.
%   - Data     : The specific data for each problem.
%   - x0       : The initial points for each test problem.
%   - alpha0   : A positive double.
%   - beta     : A double in the range (0,1).
%   - rho      : A double in the range (0,1).
%
% Outputs:
%   - H        : An array of doubles representing the minimal of the found function values for each method plus SA on each test problem.

    MDGP = 0;
    [N,R] = size(x0);
    [~,K] = size(Methods);
    H = zeros(N*R,K+1);
    
    for i = 1:N
        for j = 1:R
            index = (i-1)*R+j;
            for k = 1:K
                [List,~] = Methods{k}(MDGP,Problems{i}{1},Data{i},x0{i,j},alpha0,beta,rho,Problems{i}{3},Param{k});
                H(index,k) = min(List);
            end
            options = optimoptions(@simulannealbnd,'MaxIterations',10^10,'ObjectiveLimit',10^(-16),'FunctionTolerance',10^(-16),'MaxFunctionEvaluations',100*(Problems{i}{2}+1));
            x = x0{i,j};
            fun = @(x)Problems{i}{1}("function",x,Data{i});
            [~,fval] = simulannealbnd(fun,x,[],[],options);
            H(index,K+1) = fval;
        end
    end
end

function H = SA_MDGP(Methods,Param,Problems,Data,x0,alpha0,beta,rho)
% The function SA_MDGP
%
% Parameters:
%   - Methods  : The desired methods to compare.
%   - Param    : The specific parameters for each method.
%   - Problems : The MDGP problems.
%   - Data     : The specific data for MDGP.
%   - x0       : The initial points for each test problem.
%   - alpha0   : A positive double.
%   - beta     : A double in the range (0,1).
%   - rho      : A double in the range (0,1).
%
% Outputs:
%   - H        : An array of doubles representing the minimal of the found function values for each method plus SA on the MDGP.

    MDGP = 1;
    [M,~] = size(Data{1}{1});
    [N,R] = size(x0);
    [~,np] = size(x0{1,1});
    [~,K] = size(Methods);
    
    epsilon = 10^(-5);
    alpha0 = 1/M;
    
    fmin = M*epsilon;
    H = zeros(N*R,K+1);
    
    for i = 1:N
        for j = 1:R
            index = (i-1)*R+j;
            for k = 1:K
                [List,~] = Methods{k}(MDGP,Problems{i},Data{i},x0{i,j},alpha0,beta,rho,fmin,Param{k});
                H(index,k) = min(List);
            end
            options = optimoptions(@simulannealbnd,'MaxIterations',10^10,'ObjectiveLimit',10^(-16),'FunctionTolerance',10^(-16),'MaxFunctionEvaluations',100*(np+1));
            x = x0{i,j};
            fun = @(x)Problems{i}{1}("function",x,Data{i});
            [~,fval] = simulannealbnd(fun,x,[],[],options);
            H(index,K+1) = fval;
        end
    end
end

function H = Profiles(Methods,Param,Problems,Data,x0,alpha0,beta,rho)
% The function Profiles
%
% Parameters:
%   - Methods  : The desired methods to compare.
%   - Param    : The specific parameters for each method.
%   - Problems : The considered problems.
%   - Data     : The specific data for each problem.
%   - x0       : The initial points for each test problem.
%   - alpha0   : A positive double.
%   - beta     : A double in the range (0,1).
%   - rho      : A double in the range (0,1).
%
% Outputs:
%   - H        : An array of doubles representing the historic of the found function values for each method on each test problem.
    MDGP = 0;
    [N,R] = size(x0);
    [~,K] = size(Methods);
    
    max_np = 0;
    for i = 1:N
        if Problems{i}{2} > max_np
            max_np = Problems{i}{2};
        end
    end
    H = zeros(100*(max_np+1),N*R,K);
    for i = 1:N
        for j = 1:R
            index = (i-1)*R+j;
            for k = 1:K
                [H(:,index,k),~] = Methods{k}(MDGP,Problems{i}{1},Data{i},x0{i,j},alpha0,beta,rho,Problems{i}{3},Param{k});
            end
        end
    end
end

function H = Profiles_MDGP(Methods,Param,Problems,Data,x0,alpha0,beta,rho)
% The function Profiles_MDGP
%
% Parameters:
%   - Methods  : The desired methods to compare.
%   - Param    : The specific parameters for each method.
%   - Problems : The MDGP problems.
%   - Data     : The specific data for MDGP.
%   - x0       : The initial points for each test problem.
%   - alpha0   : A positive double.
%   - beta     : A double in the range (0,1).
%   - rho      : A double in the range (0,1).
%
% Outputs:
%   - H        : An array of doubles representing the historic of the found function values for each method on the MDGP.

    MDGP = 1;
    [M,~] = size(Data{1}{1});
    [N,R] = size(x0);
    [~,np] = size(x0{1,1});
    [~,K] = size(Methods);
  
    epsilon = 10^(-5);
    alpha0 = 1/M;
    
    fmin = M*epsilon;
    H = zeros(100*(np+1),N*R,K);
    succeed = zeros(N*R,K);
    for i = 1:N
        for j = 1:R
            index = (i-1)*R+j;
            for k = 1:K
                [List,k_ite] = Methods{k}(MDGP,Problems{1},Data{i},x0{i,j},alpha0,beta,rho,fmin,Param{k});
                if min(List) <= fmin
                    succeed(index,k) = 1;
                end
                H(:,index,k) = List;
                H(k_ite+1:100*(np+1),index,1) = List(end)*ones(100*(np+1)-k_ite,1);
            end
        end
    end
end







