function [G,n,m] = ErdosRenyiRandomgraph_dir(n, p, X, Y, showPlot)

%%***********************************************************************%
%*                        Erdos-Renyi random graph                      *%
%*           Generates a Erdos-Renyi random graph by probability        *%
%*                                                                      *%
%*                                                                      *%
%* Author: Preetham Manjunatha                                          *%
%* Github link: https://github.com/preethamam                           *%
%* Date: 05/14/2022                                                     *%
%************************************************************************%
%
%************************************************************************%
%
% Usage: [G,n,m]                = ErdosRenyiRandomgraph(n,p)
%  
% Inputs:
%
%           n                   - Graph size, number of vertexes, |V|
%           p                   - Probability, p, to connect vertices of Erdos-Renyi model      
%
%
% 
% Outputs: 
%
%           G                   - Graph output
%           n                   - Number of vertexes, |V|
%           m                   - Number of edges, |E|
%
%
%--------------------------------------------------------------------------

%------------------------------------------------------------------------------------------------------------------------
% nargin check
if nargin < 4
    error('Not enough input arguments.');
elseif nargin > 5
    error('Too many input arguments.');
end

%------------------------------------------------------------------------------------------------------------------------
% Generate random graph
% G = spones(triu(sprand(n,n,p),1));

% Initialize adjacency matrix
adjMatrix = zeros(n);


% Fill the adjacency matrix
for i = 1:n
    for j = 1:n
        if i ~= j
            adjMatrix(i, j) = rand() < p;
        end
    end
end




if nargout>2
    m = nnz(adjMatrix);
end

%------------------------------------------------------------------------------------------------------------------------
% Generate random graph with symmetry
% G = G + G';

G = digraph(adjMatrix);
G.Nodes.XData = X';
G.Nodes.YData = Y';

% Plot the graph
if (showPlot)    
    figure;
    plot_handle = plot(G,'XData',X,'YData',Y);
    plot_handle.MarkerSize = 3;
    title('Erdos-Renyi Graph')
    xlabel('X'); ylabel('Y');

    % Save the figure
    figureName = sprintf('ErdosRenyiGraph_p_%.2f.eps', p);
    saveas(plot_handle, figureName, 'epsc');

end


end