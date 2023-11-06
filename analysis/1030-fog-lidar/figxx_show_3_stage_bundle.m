clc;
close all;
global LISMO_VARS;

%% Parameter definition
x = [-0.06, 0.1; 0.9, 1.12; 1.92, 2.2];
y = - [0, 0; 2.45, 2.45; 4.98, 4.98];

%% Display
figure('Position', [0, 0, 300, 350], 'Units', 'Pixels', 'Color', 'w');

hold on;
plot([0, 0], [-100, 100], '--k');
plot([-100, 100], [0, 0], '--k');
p1 = plot([x(1, 1), x(1, 2)], [y(1, 1), y(1, 2)], 'LineWidth', 2, 'Color', [255, 0, 255] / 255, 'DisplayName', '5000 m');
p2 = plot([x(2, 1), x(2, 2)], [y(2, 1), y(2, 2)], 'LineWidth', 2, 'Color', [54, 116, 255] / 255, 'DisplayName', '100 m');
p3 = plot([x(3, 1), x(3, 2)], [y(3, 1), y(3, 2)], 'LineWidth', 2, 'Color', [29, 157, 228] / 255, 'DisplayName', '50 m');
hold off;

xlim([-0.5, 3.5]);
ylim([-6.5, 1]);

xlabel('x (mm)');
ylabel('y (mm)');
title(sprintf('Imagings of laser beam\nfrom different distances'));

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'Ticklen', [0.02, 0.02], 'fontsize', 11);

legend([p1, p2, p3], 'Location', 'NorthEast');

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', '3-stage-sim.png'), '-r300');