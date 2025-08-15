% Sample entropy on a Logistic/Chaotic map plot (annotated)
clear all
load("Log2.mat", "res")

% Parameter grid for bifurcation parameter rho
a1 = 3.55;
a2 = 4.00;
def = 10e-5;            % Step size for rho: 0.0001

%% Plot
figure(1)

% Plot three curves from res over the parameter grid a1:def:a2
% Enlarged classical SampEn curve with colour alterations for better
% display
h1 = plot(a1:def:a2, res(:,1), '-', 'LineWidth', 1.5, 'Color', [1 0 0]);  % Classical SampEn (red)
hold on
h2 = plot(a1:def:a2, res(:,3), '--', 'LineWidth', 1, 'Color', [0 0.4470 0.7410]); % Directed SampEn_G (blue, dashed)
hold on
plot(a1:def:a2, res(:,2), '-', 'LineWidth', 1, 'Color', [0 0 0]); % Undirected SampEn_G (black)

hold off
axis tight
set(gcf,'color','white')
set(gca,'FontSize',19)
box off

h_legend = legend('Classical SampEn','Directed SampEn_G', 'Undirected SampEn_G','Location','northwest');
xlim([0.99*3.55 4*1.01])
ylim([0 .65])
xlabel('Birfurcation parameter \rho')
ylabel('Sample Entropy value')

set(h_legend, 'FontSize', 14, 'NumColumns', 1, 'Box', 'off');

% Figure size 
fig_width = 12.6; % inches
fig_height = 3;   % inches
set(gcf, 'Units', 'Inches', 'Position', [1, 1, fig_width, fig_height], ...
    'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height], ...
    'PaperPositionMode', 'auto');

% Save figure
saveas(gca,'LogisticMap3_plot.eps','epsc');
