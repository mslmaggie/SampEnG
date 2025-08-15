clear all

%% ------------------------------------------------------------------------
% Parameters for Plot
r_list  = 0.05:0.05:0.40;   % tolerance values



load("SEG_Brodatz_texture_m2_25patch.mat")
%% ------------------------------------------------------------------------
% Plot
figure; hold on
mk = {'-ko','-rv','-gs','-yd','-b^','-mh','-co','-kv','-rs'};  % markers
for k = 1:Nimg
    plot(r_list, res(k,:), mk{k}, 'LineWidth',1.3);
end
xlabel('r')
ylabel('Sample Entropy_G value')
h_legend = legend({'05','15','30','36','45','75','93','95','102'}, ...
       'Location','northeast');

hold off

set(gca, 'FontSize', 19);
set(gcf, 'Color', 'white');

xlim([0.05,0.4])

set(h_legend, 'FontSize', 19, 'NumColumns', 3, 'Box', 'off')

fig_width = 6.3;  
fig_height = 4.2; 

set(gcf, 'Units', 'Inches', 'Position', [1, 1, fig_width, fig_height], ...
         'PaperUnits', 'Inches', 'PaperSize', [fig_width, fig_height], ...
         'PaperPositionMode', 'auto');
set(gca, 'LooseInset', get(gca, 'TightInset')); % Reduce whitespace

%% ------------------------------------------------------------------------

saveas(gcf,'SEG_Brodatz_texture_m2_25patch.eps','epsc')
