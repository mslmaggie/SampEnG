%//%************************************************************************%
%//%*                    Random graphs and trees                           *%
%//%*                                                                      *%
%//%*          Authors: Dr. Preetham Manjunatha & Zhenghao Li              *%
%//%*             GitHub: https://github.com/preethamam                    *%
%//%*                                                                      *%
%//%************************************************************************%
%//%*                                                                      *%                             
%//%*             University of Southern california,                       *%
%//%*             Los Angeles, California.                                 *%
%//%************************************************************************%

%% Start parameters
%--------------------------------------------------------------------------
clear; close all; clc;

%% Erodos Renyi
showPlotER = 0;         % Flag to display Erdos-Renyi plots

%% Pre-processing
%--------------------------------------------------------------------------
% Boundary (rectangular or square)
% imX1, imX2 define the X-axis bounds; imY1, imY2 define the Y-axis bounds
imX1 = 1; imX2 = 10;
imY1 = 1; imY2 = 10;

%--------------------------------------------------------------------------
% Define X and Y points for random graph
% (X and Y are used to randomize node positions for graph generation.)

%% Random graph generation methods

%-------------
% List of edge connectivity values to test (p in ER model)
p_list = [ .02, .03, .04, .05, .06, .08, .1, .2, .3, .4, .5];
p_length = size(p_list,2);       % Number of p values

repetitions = 20;                           % Independent trials per p value
num_m = 3;                                  % Embedding dimension m for SampEn
res = zeros(p_length, repetitions, 3);      % Results: [p × trials × m]
res_p = zeros(p_length, repetitions, 3);    % Placeholder (not used below)
t = 1;                                      % Trial counter
numPoints = 100;                            % Number of nodes in the ER graph

for j = 1:repetitions
    for i = 1:size(p_list,2)
        
        r = 0.2;  % Tolerance for SampEn matching

        % Randomized coordinates for nodes (scaled by 1/100)
        X = randi([imX1,imX2], 1, numPoints) / 100;
        Y = randi([imY1,imY2], 1, numPoints) / 100;

        p = p_list(i);  % Current edge connectivity

        % Generate a directed random graph using a custom ER function
        G1 = ErdosRenyiRandomgraph_dir(numPoints, p, X, Y, showPlotER);
        
        % Convert graph edges to an adjacency matrix (Adj)
        Adj = edgesToAdjMatrix(G1.Edges, numPoints);  % Adjacency matrix

        for m = 1:num_m
            % Set signal as a random signal Y
            signal = Y;
            
            % Compute SampEnG
            SampEn = SEG(signal, Adj, m, r, 'chebychev');

            % Store the result for this (p, trial, m)
            res(i, t, m) = SampEn;
        end
    end
    t = t + 1;  % Move to the next trial

    disp({i,j});  % Progress indicator
end



%% Averaging
% Parameters
num_p = length(p_list);        % Number of p values tested
num_trials = size(res, 2);     % Number of repetitions (trials)
m_vals = 1:num_m;               % Range of m values

% Initialize containers for statistics
mean_vals = zeros(num_p, num_m);  % Mean SampEn over trials for each p,m
std_vals  = zeros(num_p, num_m);  % Std dev over trials for each p,m

% Compute mean and std dev over repetitions for each p and m
for i = 1:num_p
    for m = 1:num_m
        sampens = squeeze(res(i, :, m));   % Extract SampEn values for this (p,m)
        mean_vals(i, m) = mean(sampens);
        std_vals(i, m)  = std(sampens);
    end
end


%% Colour map
N = 16;                      % Number of colors
cmap = hsv(N);               % Original HSV colormap
hsv_map = rgb2hsv(cmap);     % Convert to HSV
hsv_map(:,2) = hsv_map(:,2) * 0.8;  % Slightly desaturate
hsv_map(:,3) = hsv_map(:,3) * 0.79; % Slightly darken
cmap_dimmer = hsv2rgb(hsv_map);    % Back to RGB

colors = colormap(cmap_dimmer);     % Apply the dimmed colormap

%% Plotting

figure;
hold on;
for i = 1:num_p
    % Plot mean SampEn with error bars (std) across trials for each p
    errorbar(m_vals, mean_vals(i, :), std_vals(i, :), '-o', ...
        'LineWidth', 1, 'DisplayName', sprintf('p = %.2f', p_list(i)), ...
        'Color', colors(i, :));   % Color assigned from the custom colormap
end
xlabel('Embedding Dimension (m)');
ylabel('Sample Entropy_G value');  % Note: "G" denotes a variant; adjust if needed

h_legend = legend('Location', 'northeast');
set(h_legend, 'FontSize', 15, 'NumColumns', 3, 'Box', 'off');
set(gca, 'FontSize', 19);
set(gcf, 'Color', 'white');

% Figure sizing for publication
fig_width = 6.3; 
fig_height = 4.2; 

set(gcf, 'Units', 'Inches', 'Position', [1, 1, fig_width, fig_height], ...
         'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height], ...
         'PaperPositionMode', 'auto');
set(gca, 'LooseInset', get(gca, 'TightInset')); % Reduce whitespace

xticks([1 2 3]);               % Tick positions corresponding to m = 1..3
xticklabels({'1', '2', '3'});  % Tick labels

ylim([-.3,3.2])                 % Y-axis limits

% grid on;
hold off;

% Save outputs
saveas(gcf, 'SampleEntropy_ER_graph.eps', 'epsc');
save('ER_SEG_mean.mat', 'mean_vals');
save('ER_SEG_std.mat', 'std_vals');



%% Convert to Adjacency matrix format

function adjMatrix = edgesToAdjMatrix(edges, numNodes)
    % edges: a N x 3 matrix where each row is [startNode, endNode, weight]
    % numNodes: the total number of nodes in the graph

    % Initialize the adjacency matrix with zeros
    adjMatrix = zeros(numNodes);

    % Iterate over the edges and fill the adjacency matrix
    for i = 1:size(edges, 1)
        startNode = edges.EndNodes(i, 1);
        endNode   = edges.EndNodes(i, 2);
        weight    = edges.Weight(i,1);
        adjMatrix(startNode, endNode) = weight;
    end
end
