clc;
close all;
projectDir = fileparts(fileparts(mfilename('fullpath')));

%% Parameter definition
kernelFile = fullfile(projectDir, 'data', 'mie_kernel_water_droplets.mat');

%% Mie scattering
data = load(kernelFile, 'wavelength', 'r', 'qe', 'qs', 'qb', 'g');
qe = data.qe;
qs = data.qs;
qb = data.qb;
g = data.g;
r = data.r;
wavelength = data.wavelength;

%% Fog size dist.
fogSDAdvMod = fogSD_SS('type', 'advection', 'intensity', 'moderate');
fogSDAdvHvy = fogSD_SS('type', 'advection', 'intensity', 'heavy');
fogSDRadMod = fogSD_SS('type', 'radiation', 'intensity', 'moderate');
fogSDRadHvy = fogSD_SS('type', 'radiation', 'intensity', 'heavy');

backscatter = NaN(4, length(wavelength));
extinction = NaN(4, length(wavelength));
asymFactor = NaN(4, length(wavelength));
for iWL = 1:length(wavelength)
    backscatter(1, iWL) = sum(fogSDAdvMod(r) .* qb(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1))) / (4 * pi);
    extinction(1, iWL) = sum(fogSDAdvMod(r) .* qe(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1)));
    asymFactor(1, iWL) = sum(fogSDAdvMod(r) .* qs(iWL, :) .* g(iWL, :) .* r.^2 * (r(2) - r(1))) ./ sum(fogSDAdvMod(r) .* qs(iWL, :) .* r.^2 * (r(2) - r(1)));

    backscatter(2, iWL) = sum(fogSDAdvHvy(r) .* qb(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1))) / (4 * pi);
    extinction(2, iWL) = sum(fogSDAdvHvy(r) .* qe(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1)));
    asymFactor(2, iWL) = sum(fogSDAdvHvy(r) .* qs(iWL, :) .* g(iWL, :) .* r.^2 * (r(2) - r(1))) ./ sum(fogSDAdvHvy(r) .* qs(iWL, :) .* r.^2 * (r(2) - r(1)));

    backscatter(3, iWL) = sum(fogSDRadMod(r) .* qb(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1))) / (4 * pi);
    extinction(3, iWL) = sum(fogSDRadMod(r) .* qe(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1)));
    asymFactor(3, iWL) = sum(fogSDRadMod(r) .* qs(iWL, :) .* g(iWL, :) .* r.^2 * (r(2) - r(1))) ./ sum(fogSDRadMod(r) .* qs(iWL, :) .* r.^2 * (r(2) - r(1)));

    backscatter(4, iWL) = sum(fogSDRadHvy(r) .* qb(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1))) / (4 * pi);
    extinction(4, iWL) = sum(fogSDRadHvy(r) .* qe(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1)));
    asymFactor(4, iWL) = sum(fogSDRadHvy(r) .* qs(iWL, :) .* g(iWL, :) .* r.^2 * (r(2) - r(1))) ./ sum(fogSDRadHvy(r) .* qs(iWL, :) .* r.^2 * (r(2) - r(1)));
end

%% Output results
fprintf('Extinction coefficient of clouds at different wavelength (km-1)...\n');
factor = 6.5 / interp1(wavelength, extinction(2, :), 0.355);
fprintf('@355nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, extinction(2, :) * factor, 0.355), 4 / (interp1(wavelength, extinction(2, :) * factor, 0.355) + 0.0124));
fprintf('@532nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, extinction(2, :) * factor, 0.532), 4 / (interp1(wavelength, extinction(2, :) * factor, 0.532) + 0.0024));
fprintf('@1064nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, extinction(2, :) * factor, 1.064), 4 / interp1(wavelength, extinction(2, :) * factor, 1.064));
fprintf('@1600nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, extinction(2, :) * factor, 1.6), 4 / interp1(wavelength, extinction(2, :) * factor, 1.6));
fprintf('@2000nm: %5.3f; penatration depth (km): %5.3f\n', interp1(wavelength, extinction(2, :) * factor, 2.0), 4 / interp1(wavelength, extinction(2, :) * factor, 2.0));

%% Data visualization
subPos = subfigPos([0.12, 0.10, 0.85, 0.84], 3, 1, 0, 0.04);

figure('Position', [0, 10, 400, 500], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', subPos(1, :), 'Units', 'normalized');
p2 = plot(wavelength, extinction(2, :), 'LineStyle', '-', 'Color', [211, 211, 211]/255, 'LineWidth', 2, 'DisplayName', 'Severe Advection fog'); hold on;
p3 = plot(wavelength, extinction(3, :), 'LineStyle', '--', 'Color', [65, 105, 226]/255, 'LineWidth', 2, 'DisplayName', 'Moderate Radiation fog'); hold on;
p4 = plot(wavelength, extinction(4, :), 'LineStyle', '--', 'Color', [135, 206, 250]/255, 'LineWidth', 2, 'DisplayName', 'Severe Radiation fog'); hold on;
p1 = plot(wavelength, extinction(1, :), 'LineStyle', '-', 'Color', [0, 0, 0]/255, 'LineWidth', 2, 'DisplayName', 'Moderate Advection fog'); hold on;

xlabel('');
ylabel('Extinction (km^{-1})');
title('Optical properties for moderate advection fog');

xlim([wavelength(1), wavelength(end)]);
ylim([0, 50]);


set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'XTickLabel', '', 'Box', 'on', 'TickLen', [0.02, 0.01]);

subplot('Position', subPos(2, :), 'Units', 'normalized');
p2 = plot(wavelength, backscatter(2, :), 'LineStyle', '-', 'Color', [211, 211, 211]/255, 'LineWidth', 2, 'DisplayName', 'Severe Advection fog'); hold on;
p3 = plot(wavelength, backscatter(3, :), 'LineStyle', '--', 'Color', [65, 105, 226]/255, 'LineWidth', 2, 'DisplayName', 'Moderate Radiation fog'); hold on;
p4 = plot(wavelength, backscatter(4, :), 'LineStyle', '--', 'Color', [135, 206, 250]/255, 'LineWidth', 2, 'DisplayName', 'Severe Radiation fog'); hold on;
p1 = plot(wavelength, backscatter(1, :), 'LineStyle', '-', 'Color', [0, 0, 0]/255, 'LineWidth', 2, 'DisplayName', 'Moderate Advection fog'); hold on;

xlabel('');
ylabel('Backscatter (km^{-1}sr^{-1})');

xlim([wavelength(1), 3]);
ylim([0, 3]);

legend([p1, p2, p3, p4], 'Location', 'NorthEast');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'XTickLabel', '', 'Box', 'on', 'TickLen', [0.02, 0.01]);

subplot('Position', subPos(3, :), 'Units', 'normalized');
% p2 = plot(wavelength, extinction(2, :) ./ backscatter(2, :), 'LineStyle', '-', 'Color', [211, 211, 211]/255, 'LineWidth', 2, 'DisplayName', 'Severe Advection fog'); hold on;
% p3 = plot(wavelength, extinction(3, :) ./ backscatter(3, :), 'LineStyle', '--', 'Color', [65, 105, 226]/255, 'LineWidth', 2, 'DisplayName', 'Moderate Radiation fog'); hold on;
% p4 = plot(wavelength, extinction(4, :) ./ backscatter(4, :), 'LineStyle', '--', 'Color', [135, 206, 250]/255, 'LineWidth', 2, 'DisplayName', 'Severe Radiation fog'); hold on;
p1 = plot(wavelength, extinction(1, :) ./ backscatter(1, :), 'LineStyle', '-', 'Color', [0, 0, 0]/255, 'LineWidth', 2, 'DisplayName', 'Moderate Advection fog'); hold on;

xlabel('\lambda (\mum)');
ylabel('Lidar Ratio (sr)');

xlim([wavelength(1), 3]);
ylim([0, 50]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.01]);

export_fig(gcf, fullfile(projectDir, 'image', 'ensemble_sea_fog_optical_properties.png'), '-r300');

lrAdvModerate = extinction(1, :) ./ backscatter(1, :);
lrAdvHeavy = extinction(2, :) ./ backscatter(2, :);
lrRadModerate = extinction(3, :) ./ backscatter(3, :);
lrRadHeavy = extinction(4, :) ./ backscatter(4, :);
save(fullfile(projectDir, 'data', 'sea-fog-lidar-ratio.mat'), 'wavelength', 'lrAdvModerate', 'lrAdvHeavy', 'lrRadModerate', 'lrRadHeavy');

%% asymmetry factor
figure('Position', [0, 10, 400, 300], 'Units', 'Pixels', 'Color', 'w');

p2 = plot(wavelength, asymFactor(2, :), 'LineStyle', '-', 'Color', [211, 211, 211]/255, 'LineWidth', 2, 'DisplayName', 'Severe Advection fog'); hold on;
p3 = plot(wavelength, asymFactor(3, :), 'LineStyle', '--', 'Color', [65, 105, 226]/255, 'LineWidth', 2, 'DisplayName', 'Moderate Radiation fog'); hold on;
p4 = plot(wavelength, asymFactor(4, :), 'LineStyle', '--', 'Color', [135, 206, 250]/255, 'LineWidth', 2, 'DisplayName', 'Severe Radiation fog'); hold on;
p1 = plot(wavelength, asymFactor(1, :), 'LineStyle', '-', 'Color', [0, 0, 0]/255, 'LineWidth', 2, 'DisplayName', 'Moderate Advection fog'); hold on;

xlabel('\lambda (\mum)');
ylabel('g');

xlim([wavelength(1), wavelength(end)]);
ylim([0, 1]);
legend([p1, p2, p3, p4], 'Location', 'SouthEast');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.01]);

export_fig(gcf, fullfile(projectDir, 'image', 'ensemble_sea_fog_asym_factor.png'), '-r300');