clc;
close all;

%% Parameter Definition
kernerlFile = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data', 'mie_kernel_water_droplets.mat');

%% Load Mie Kernels
data = load(kernerlFile, 'wavelength', 'r', 'qe', 'qs', 'qb', 'g');
qe = data.qe;
qs = data.qs;
qb = data.qb;
g = data.g;
r = data.r;
wavelength = data.wavelength;

%% Cloud Droplet Size Distribution
marineSLCD = slwcDSD_mgd('N0', 7.4e6, 'gamma0', 8.6, 'Rm', 1.35, 'type', 'input');
continentalSLCD = slwcDSD_mgd('N0', 28.8e6, 'gamma0', 8.7, 'Rm', 0.65, 'type', 'input');

%% Cloud Optical Properties
marineSLBackscatter = NaN(length(wavelength), 1);
marineSLExtinction = NaN(length(wavelength), 1);
marineSLAsymFactor = NaN(length(wavelength), 1);
continentalSLBackscatter = NaN(length(wavelength), 1);
continentalSLExtinction = NaN(length(wavelength), 1);
continentalSLAsymFactor = NaN(length(wavelength), 1);

for iWL = 1:length(wavelength)
    marineSLBackscatter(iWL) = sum(marineSLCD(r) .* qb(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1))) / (4 * pi);
    marineSLExtinction(iWL) = sum(marineSLCD(r) .* qe(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1)));
    marineSLAsymFactor(iWL) = sum(marineSLCD(r) .* qs(iWL, :) .* g(iWL, :) .* r.^2 * (r(2) - r(1))) ./ sum(marineSLCD(r) .* qs(iWL, :) .* r.^2 * (r(2) - r(1)));

    continentalSLBackscatter(iWL) = sum(continentalSLCD(r) .* qb(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1))) / (4 * pi);
    continentalSLExtinction(iWL) = sum(continentalSLCD(r) .* qe(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1)));
    continentalSLAsymFactor(iWL) = sum(continentalSLCD(r) .* qs(iWL, :) .* g(iWL, :) .* r.^2 * (r(2) - r(1))) ./ sum(continentalSLCD(r) .* qs(iWL, :) .* r.^2 * (r(2) - r(1)));
end

%% Output Results
fprintf('Extinction coefficient of marine stratiform low-level clouds at different wavelength (km-1)...\n');
factor = 1;
fprintf('@355nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, marineSLExtinction * factor, 0.355), 4 / (interp1(wavelength, marineSLExtinction * factor, 0.355) + 0.0124));
fprintf('@532nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, marineSLExtinction * factor, 0.532), 4 / (interp1(wavelength, marineSLExtinction * factor, 0.532) + 0.0024));
fprintf('@1064nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, marineSLExtinction * factor, 1.064), 4 / interp1(wavelength, marineSLExtinction * factor, 1.064));
fprintf('@1600nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, marineSLExtinction * factor, 1.6), 4 / interp1(wavelength, marineSLExtinction * factor, 1.6));
fprintf('@2000nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, marineSLExtinction * factor, 2.0), 4 / interp1(wavelength, marineSLExtinction * factor, 2.0));

fprintf('\nExtinction coefficient of continental stratiform low-level clouds at different wavelength (km-1)...\n');
factor = 1;
fprintf('@355nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, continentalSLExtinction * factor, 0.355), 4 / (interp1(wavelength, continentalSLExtinction * factor, 0.355) + 0.0124));
fprintf('@532nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, continentalSLExtinction * factor, 0.532), 4 / (interp1(wavelength, continentalSLExtinction * factor, 0.532) + 0.0024));
fprintf('@1064nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, continentalSLExtinction * factor, 1.064), 4 / interp1(wavelength, continentalSLExtinction * factor, 1.064));
fprintf('@1600nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, continentalSLExtinction * factor, 1.6), 4 / interp1(wavelength, continentalSLExtinction * factor, 1.6));
fprintf('@2000nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, continentalSLExtinction * factor, 2.0), 4 / interp1(wavelength, continentalSLExtinction * factor, 2.0));

%% Data Visualization

% microphysics for stratiform low-level clouds
figure('Position', [0, 10, 500, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(r, marineSLCD(r) * 1e-6, 'LineStyle', '-', 'Color', [0, 0, 0]/255, 'LineWidth', 2, 'DisplayName', 'Marine Clouds');
p2 = plot(r, continentalSLCD(r) * 1e-6, 'LineStyle', '-', 'Color', [211, 211, 211]/255, 'LineWidth', 2, 'DisplayName', 'Continental Clouds');
hold off;

xlabel('Radius ( \mum )');
ylabel('N (cm^{-3} \mum^{-1})');

xlim([0, 60]);
ylim([0, 10]);
legend([p1, p2], 'Location', 'NorthEast');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.01]);

export_fig(gcf, fullfile(projectDir, 'image', 'ensemble_low-level-clouds_microphysics.png'), '-r300');

% optical properties for stratiform low-level clouds
figure('Position', [0, 10, 400, 500], 'Units', 'Pixels', 'Color', 'w');

subPos = subfigPos([0.15, 0.10, 0.78, 0.84], 3, 1, 0, 0.04);

subplot('Position', subPos(1, :), 'Units', 'normalized');

hold on;
plot(wavelength, continentalSLExtinction, 'LineStyle', '-', 'Color', [211, 211, 211]/255, 'LineWidth', 2, 'DisplayName', 'Continental Clouds');
plot(wavelength, marineSLExtinction, 'LineStyle', '-', 'Color', [0, 0, 0]/255, 'LineWidth', 2, 'DisplayName', 'Marine Clouds');
hold off;

xlabel('');
ylabel('Extinction (km^{-1})');
title('Optical properties for stratiform low-level clouds');

xlim([wavelength(1), wavelength(end)]);
ylim([0, 10]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'XTickLabel', '', 'Box', 'on', 'TickLen', [0.02, 0.01], 'YScale', 'Linear', 'FontSize', 11);

subplot('Position', subPos(2, :), 'Units', 'normalized');

hold on;
plot(wavelength, continentalSLBackscatter, 'LineStyle', '-', 'Color', [211, 211, 211]/255, 'LineWidth', 2, 'DisplayName', 'Continental Clouds');
plot(wavelength, marineSLBackscatter, 'LineStyle', '-', 'Color', [0, 0, 0]/255, 'LineWidth', 2, 'DisplayName', 'Marine Clouds');
hold off;

xlabel('');
ylabel('Backscatter (km^{-1}sr^{-1})');

xlim([wavelength(1), wavelength(end)]);
ylim([1e-3, 1e1]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'XTickLabel', '', 'Box', 'on', 'TickLen', [0.02, 0.01], 'YScale', 'log', 'YTick', 10.^(-3:1), 'FontSize', 11);

subplot('Position', subPos(3, :), 'Units', 'normalized');

hold on;
p2 = plot(wavelength, continentalSLExtinction ./ continentalSLBackscatter, 'LineStyle', '-', 'Color', [211, 211, 211]/255, 'LineWidth', 2, 'DisplayName', 'Continental Clouds');
p1 = plot(wavelength, marineSLExtinction ./ marineSLBackscatter, 'LineStyle', '-', 'Color', [0, 0, 0]/255, 'LineWidth', 2, 'DisplayName', 'Marine Clouds');
hold off;

xlabel('\lambda (\mum)');
ylabel('Lidar Ratio (sr)');

xlim([wavelength(1), wavelength(end)]);
ylim([1e0, 1e4]);

legend([p1, p2], 'Location', 'NorthWest');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.01], 'YScale', 'log', 'FontSize', 11, 'YTick', [1, 10, 100], 'YTickLabel', {'1', '10', '100'});

export_fig(gcf, fullfile(projectDir, 'image', 'ensemble_low-level-clouds_optical_properties.png'), '-r300');