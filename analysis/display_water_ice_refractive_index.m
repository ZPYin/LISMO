clc;
close all;
projectDir = fileparts(fileparts(mfilename('fullpath')));

%% Parameters Definition
wopFile = 'refractive_index_water.mat';
iopFile = 'IOP_2008_ASCIItable.dat';

%% Load Data

% Load refractive index of water
wop = load(fullfile(projectDir, 'data', wopFile));

% Load inherent optical properties of water
fid = fopen(fullfile(projectDir, 'data', iopFile), 'r');

data = textscan(fid, '%f%f%f', 'HeaderLines', 4, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);

fclose(fid);

%% Output results
fprintf('Complex refractive index of water...\n');
fprintf('@355nm: %5.3f - %4.2e i\n', interp1(wop.wavelength, wop.refractive_index_real, 0.355), interp1(wop.wavelength, wop.refractive_index_imaginary, 0.355));
fprintf('@532nm: %5.3f - %4.2e i\n', interp1(wop.wavelength, wop.refractive_index_real, 0.532), interp1(wop.wavelength, wop.refractive_index_imaginary, 0.532));
fprintf('@1064nm: %5.3f - %4.2e i\n', interp1(wop.wavelength, wop.refractive_index_real, 1.064), interp1(wop.wavelength, wop.refractive_index_imaginary, 1.064));
fprintf('@1600nm: %5.3f - %4.2e i\n', interp1(wop.wavelength, wop.refractive_index_real, 1.6), interp1(wop.wavelength, wop.refractive_index_imaginary, 1.6));
fprintf('@2000nm: %5.3f - %4.2e i\n', interp1(wop.wavelength, wop.refractive_index_real, 2.0), interp1(wop.wavelength, wop.refractive_index_imaginary, 2.0));

fprintf('Complex refractive index of ice...\n');
fprintf('@355nm: %5.3f - %4.2e i\n', interp1(data{1}, data{2}, 0.355), interp1(data{1}, data{3}, 0.355));
fprintf('@532nm: %5.3f - %4.2e i\n', interp1(data{1}, data{2}, 0.532), interp1(data{1}, data{3}, 0.532));
fprintf('@1064nm: %5.3f - %4.2e i\n', interp1(data{1}, data{2}, 1.064), interp1(data{1}, data{3}, 1.064));
fprintf('@1600nm: %5.3f - %4.2e i\n', interp1(data{1}, data{2}, 1.6), interp1(data{1}, data{3}, 1.6));
fprintf('@2000nm: %5.3f - %4.2e i\n', interp1(data{1}, data{2}, 2.0), interp1(data{1}, data{3}, 2.0));

%% Data Display

% Refractive index of water
figure('Position', [0, 0, 700, 300], 'Units', 'Pixels', 'Color', 'w');

subPos = subfigPos([0.1, 0.14, 0.88, 0.83], 1, 2, 0.10, 0);

subplot('Position', subPos(1, :), 'Units', 'Normalized');

scatter(wop.wavelength, wop.refractive_index_real, 10, 'Marker', '^', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); hold on;
plot(wop.wavelength, wop.refractive_index_real, '-b');

xlim([0, 3]);
ylim([1.1, 1.5]);

% xlabel('\lambda (\mum)');
xlabel('入射波长 (\mum)');
% ylabel('n');
ylabel('复折射率实部');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLength', [0.02, 0.01]);

subplot('Position', subPos(2, :), 'Units', 'Normalized');

semilogy(wop.wavelength, wop.refractive_index_imaginary, 'MarkerSize', 3, 'LineStyle', 'none', 'Marker', '^', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); hold on;
plot(wop.wavelength, wop.refractive_index_imaginary, '-b');

xlim([0, 3]);
ylim([1e-10, 1e0]);

% xlabel('\lambda (\mum)');
xlabel('入射波长 (\mum)');
% ylabel('k');
ylabel('复折射率虚部');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLength', [0.02, 0.01]);

export_fig(gcf, fullfile(projectDir, 'image', 'refractive_index_of_water.png'), '-r300');

% Refractive index of ice
figure('Position', [0, 0, 700, 300], 'Units', 'Pixels', 'Color', 'w');

subPos = subfigPos([0.1, 0.14, 0.88, 0.83], 1, 2, 0.10, 0);

subplot('Position', subPos(1, :), 'Units', 'Normalized');

scatter(data{1}, data{2}, 10, 'Marker', '^', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); hold on;
plot(data{1}, data{2}, '-b');

xlim([0, 3]);
ylim([1.1, 1.5]);

% xlabel('\lambda (\mum)');
xlabel('入射波长 (\mum)');
% ylabel('n');
ylabel('复折射率实部');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLength', [0.02, 0.01]);

subplot('Position', subPos(2, :), 'Units', 'Normalized');

semilogy(data{1}, data{3}, 'MarkerSize', 3, 'LineStyle', 'none', 'Marker', '^', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); hold on;
plot(data{1}, data{3}, '-b');

xlim([0, 3]);
ylim([1e-10, 1e0]);

% xlabel('\lambda (\mum)');
xlabel('入射波长 (\mum)');
% ylabel('k');
ylabel('复折射率虚部');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLength', [0.02, 0.01]);

export_fig(gcf, fullfile(projectDir, 'image', 'refractive_index_of_ice.png'), '-r300');