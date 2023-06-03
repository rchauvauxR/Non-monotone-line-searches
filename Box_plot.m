function Box_plot(H,dic,names,M)
% The function Box_plot
% Parameters:
%   - H      : The data to plot
%   - dic    : The names of the methods used
%   - names  : The names of the problems considered
%   - M      : The number of initial points consider per problems

    [np,~] = size(H); 
    N = floor(np/M);
    for i = 1:N
        boxPlot(H((i-1)*M+1:i*M,:),1,1,dic,names(i),i);
    end
end


