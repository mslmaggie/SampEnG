
%% Load readings and geometric information from dataset
clear all
close all
S = load('grasp_molene_data.mat');
D = S.molene_graphs{2,1};

% Data name and coordinates
D.data_name          % e.g., {'temperature'} {'maximal temperature'} {'minimal temperature'}
G.coords = D.coordinates;  % Coordinates x,y,z (not used directly below)
Mat2 = D.distances;          % distance_matrix (37x37)

%% Generate undirected, weighted Adjacency matrix

Temp = D.data{2,1} - 273.15;   % Temperature in Celsius (subtract 273.15 from Kelvin)
N = 37;                         % Number of vertices/nodes
M = Mat2;                       % Start with the distance matrix as a base
sigma1 = 1.0e+05;               % Threshold for distance (large distance cutoff)
sigma2 = 5.1e+09;               % Scale parameter for Gaussian-like kernel

% Build a Adjacency matrix from distances with a distance cutoff
for ii = 1:N
    for jj = 1:N
        if Mat2(ii,jj) >= sigma1
           M(ii,jj) = 0;                 % Too far apart: set affinity to 0
        else
           M(ii,jj) = exp(-(Mat2(ii,jj))^2 / (2*sigma2));  % Weight assigned from coordinates
        end  
    end 
end 

Adj = M - eye(N);  % Adjacency matrix for a graph: weights off-diagonal only; diagonal becomes 0

%% Main loop

dist_type = 'chebychev';    % Distance metric for SEG

r_list = 0.05:0.01:0.5;     % Tolerance levels
r_length = size(r_list,2);
m = 1;                      % Embedding dimension for SampEn

for i = 1:r_length
    r = r_list(i);
    for qq = 1:744
        
        x = Temp(:, qq)';  % For each hour (column), take the 37-node temperature vector
        res(qq, i) = SEG(x, Adj, m, r, dist_type);
        
    end 
end

%% Extract results from 4AM and 2PM only

% Calculate the total number of days (recordings/hour)
num_days = 744 / 24;

% Variables to store results from 4AM and 2PM
res_4am = zeros(num_days, r_length);
res_2pm = zeros(num_days, r_length);

for day = 1:num_days
    start_idx = (day - 1) * 24;
    res_4am(day,:) = res(start_idx + 4,:);   % 4 AM = 4th hour
    res_2pm(day,:) = res(start_idx + 14,:);  % 2 PM = 14th hour
end

% Compute mean and standard error
mean_4am = mean(res_4am, 1);
mean_2pm = mean(res_2pm, 1);

se_4am = std(res_4am) / sqrt(num_days);     % Standard Error
se_2pm = std(res_2pm) / sqrt(num_days);

%% Plot

figure;
hold on;

errorbar(r_list, mean_4am, se_4am, '-o', 'DisplayName', 'Temperature recording at 4 AM', 'LineWidth', 1.5);
errorbar(r_list, mean_2pm, se_2pm, '-s', 'DisplayName', 'Temperature recording at 2 PM', 'LineWidth', 1.5);

xlabel('r');
ylabel('Sample Entropy_G value');
% title('Sample Entropy VS Tolerance_{r} for Temperature measurements', 'FontSize', 20);
h_legend = legend('show');


set(h_legend, 'FontSize', 17, 'NumColumns', 1, 'Box', 'off');
% set(h_legend, 'Orientation', 'horizontal');  % optional

set(gca, 'FontSize', 19);
set(gcf, 'Color', 'white');


fig_width = 6.3;  % inches, for two-column width in LaTeX
fig_height = 4.2; % inches, adjust as needed for your data

set(gcf, 'Units', 'Inches', 'Position', [1, 1, fig_width, fig_height], ...
         'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height], ...
         'PaperPositionMode', 'auto');
set(gca, 'LooseInset', get(gca, 'TightInset')); % Reduce whitespace

xlim([0.05,0.32])
ylim([0.2,2.1])

saveas(gcf, 'SampleEntropy_Temperature_measurement.eps', 'epsc');
