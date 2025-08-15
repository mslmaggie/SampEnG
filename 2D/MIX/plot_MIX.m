load("m_2_2D_MIX_averaged.mat")

%% Plot with inverted colormap + styled lines
figure;
hold on;

q_list = 0.1:0.1:0.9;
q_length = size(q_list,2);
vector=10:10:150; %change2

%% Define colormap 
N = 12;                      % Number of colors
cmap = hsv(N);               % Original HSV colormap
hsv_map = rgb2hsv(cmap);     % Convert to HSV
hsv_map(:,2) = hsv_map(:,2) * 0.8;  % Reduce brightness (V) to 70%
hsv_map(:,3) = hsv_map(:,3) * 0.79;  % Reduce brightness (V) to 70%
cmap_dimmer = hsv2rgb(hsv_map);    % Convert back to RGB

colormap(cmap_dimmer);       % Apply muted HSV colormap



%% Plot
lineStyles = {'-', '--', ':', '-.'};  
hList = gobjects(q_length, 1);        

for jj = 1:q_length
    data = RES1{jj};
    y = mean(data);
    errors = std(data);
    x = vector;

    % Assign inverted color and cycling style
    color_j = (cmap_dimmer(jj, :));
    style_j = lineStyles{mod(jj - 1, length(lineStyles)) + 1};

    % Plot error bar
    h = errorbar(x, y, errors, 'LineWidth', 1.2, ...
                 'Color', color_j, 'LineStyle', style_j);
    hList(jj) = h;
end

%% Formatting
ylim([0.5,1.6])

xlabel('Image size');
ylabel('Sample Entropy_{G} value');

legend_entries = arrayfun(@(q) sprintf('p = %.1f', q), q_list, 'UniformOutput', false);
h_legend = legend(hList, legend_entries, 'Location', 'southeast');

set(h_legend, 'FontSize', 18, 'NumColumns', 3, 'Box', 'off');

set(gca, 'FontSize', 19);
set(gcf, 'Color', 'white');


fig_width = 6.3;  
fig_height = 4.2; 

set(gcf, 'Units', 'Inches', 'Position', [1, 1, fig_width, fig_height], ...
         'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height], ...
         'PaperPositionMode', 'auto');
set(gca, 'LooseInset', get(gca, 'TightInset')); % Reduce whitespace


%% Save outputs
saveas(gca, 'm_2_2D_MIX_output_fit.eps', 'epsc');
