clc;
close all;

load('D:\Coding\Matlab\LISMO\data\refractive_index_water.mat');

figure('Position', [0, 0, 700, 300], 'Units', 'Pixels', 'Color', 'w');

subPos = subfigPos([0.07, 0.14, 0.88, 0.83], 1, 2, 0.10, 0);

subplot('Position', subPos(1, :), 'Units', 'Normalized');

scatter(wavelength, refractive_index_real, 10, 'Marker', '^', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); hold on;
plot(wavelength, refractive_index_real, '-b');

xlim([0, 12]);
ylim([1.1, 1.5]);

xlabel('\lambda (\mum)');
ylabel('n');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLength', [0.02, 0.01]);

subplot('Position', subPos(2, :), 'Units', 'Normalized');

semilogy(wavelength, refractive_index_imaginary, 'MarkerSize', 3, 'LineStyle', 'none', 'Marker', '^', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); hold on;
plot(wavelength, refractive_index_imaginary, '-b');

xlim([0, 12]);
ylim([1e-10, 1e0]);

xlabel('\lambda (\mum)');
ylabel('k');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLength', [0.02, 0.01]);

export_fig(gcf, 'D:\Coding\Matlab\LISMO\image\refractive_index_of_water.png', '-r300');