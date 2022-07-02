clc;
close all;

%% Parameter definition
wavelength = 0.3:0.02:4;
r = 0.01:0.1:60;

%% Mie scattering
qe = zeros(length(wavelength), length(r));
qb = zeros(length(wavelength), length(r));
for iWL = 1:length(wavelength)
    fprintf('Finished %6.2f%%\n', (iWL - 1) / length(wavelength) * 100);
    [refIdxReal, refIdxImg] = refIdxWater(wavelength(iWL) * 1e3);

    for iR = 1:length(r)
        m = refIdxReal + 1i * refIdxImg;
        k0 = 2 * pi / (wavelength(iWL) * 1e3);
        a = r(iR) * 1000;

        res = Mie(m, k0 * a);

        qe(iWL, iR) = res(1);
        qb(iWL, iR) = res(4);
    end
end

%% Fog size dist.
fogSDAdvMod = fogSD_SS('type', 'advection', 'intensity', 'moderate');
fogSDAdvHvy = fogSD_SS('type', 'advection', 'intensity', 'heavy');
fogSDRadMod = fogSD_SS('type', 'radiation', 'intensity', 'moderate');
fogSDRadHvy = fogSD_SS('type', 'radiation', 'intensity', 'heavy');

backscatter = NaN(4, length(wavelength));
extinction = NaN(4, length(wavelength));
for iWL = 1:length(wavelength)
    backscatter(1, iWL) = sum(fogSDAdvMod(r) .* qb(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1))) / (4 * pi);
    extinction(1, iWL) = sum(fogSDAdvMod(r) .* qe(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1)));

    backscatter(2, iWL) = sum(fogSDAdvHvy(r) .* qb(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1))) / (4 * pi);
    extinction(2, iWL) = sum(fogSDAdvHvy(r) .* qe(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1)));

    backscatter(3, iWL) = sum(fogSDRadMod(r) .* qb(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1))) / (4 * pi);
    extinction(3, iWL) = sum(fogSDRadMod(r) .* qe(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1)));

    backscatter(4, iWL) = sum(fogSDRadHvy(r) .* qb(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1))) / (4 * pi);
    extinction(4, iWL) = sum(fogSDRadHvy(r) .* qe(iWL, :) .* pi .* r.^2 * 1e-9 * (r(2) - r(1)));
end

%% Data visualization
subPos = subfigPos([0.12, 0.10, 0.85, 0.84], 3, 1, 0, 0.04);

figure('Position', [0, 10, 400, 500], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', subPos(1, :), 'Units', 'normalized');
p2 = plot(wavelength, extinction(2, :), 'LineStyle', '-', 'Color', [211, 211, 211]/255, 'LineWidth', 2, 'DisplayName', 'Severe Advection fog'); hold on;
p3 = plot(wavelength, extinction(3, :), 'LineStyle', '--', 'Color', [65, 105, 226]/255, 'LineWidth', 2, 'DisplayName', 'Moderate Advection fog'); hold on;
p4 = plot(wavelength, extinction(4, :), 'LineStyle', '--', 'Color', [135, 206, 250]/255, 'LineWidth', 2, 'DisplayName', 'Severe Advection fog'); hold on;
p1 = plot(wavelength, extinction(1, :), 'LineStyle', '-', 'Color', [0, 0, 0]/255, 'LineWidth', 2, 'DisplayName', 'Moderate Advection fog'); hold on;

xlabel('');
ylabel('Extinction (km^{-1})');
title('Optical properties for moderate advection fog');

xlim([wavelength(1), wavelength(end)]);
ylim([0, 50]);


set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'XTickLabel', '', 'Box', 'on', 'TickLen', [0.02, 0.01]);

subplot('Position', subPos(2, :), 'Units', 'normalized');
p2 = plot(wavelength, backscatter(2, :), 'LineStyle', '-', 'Color', [211, 211, 211]/255, 'LineWidth', 2, 'DisplayName', 'Severe Advection fog'); hold on;
p3 = plot(wavelength, backscatter(3, :), 'LineStyle', '--', 'Color', [65, 105, 226]/255, 'LineWidth', 2, 'DisplayName', 'Moderate Advection fog'); hold on;
p4 = plot(wavelength, backscatter(4, :), 'LineStyle', '--', 'Color', [135, 206, 250]/255, 'LineWidth', 2, 'DisplayName', 'Severe Advection fog'); hold on;
p1 = plot(wavelength, backscatter(1, :), 'LineStyle', '-', 'Color', [0, 0, 0]/255, 'LineWidth', 2, 'DisplayName', 'Moderate Advection fog'); hold on;

xlabel('');
ylabel('Backscatter (km^{-1}sr^{-1})');

xlim([wavelength(1), wavelength(end)]);
ylim([0, 3]);

legend([p1, p2, p3, p4], 'Location', 'NorthEast');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'XTickLabel', '', 'Box', 'on', 'TickLen', [0.02, 0.01]);

subplot('Position', subPos(3, :), 'Units', 'normalized');
p2 = plot(wavelength, extinction(2, :) ./ backscatter(2, :), 'LineStyle', '-', 'Color', [211, 211, 211]/255, 'LineWidth', 2, 'DisplayName', 'Mod. Advection fog'); hold on;
p3 = plot(wavelength, extinction(3, :) ./ backscatter(3, :), 'LineStyle', '--', 'Color', [65, 105, 226]/255, 'LineWidth', 2, 'DisplayName', 'Mod. Advection fog'); hold on;
p4 = plot(wavelength, extinction(4, :) ./ backscatter(4, :), 'LineStyle', '--', 'Color', [135, 206, 250]/255, 'LineWidth', 2, 'DisplayName', 'Mod. Advection fog'); hold on;
p1 = plot(wavelength, extinction(1, :) ./ backscatter(1, :), 'LineStyle', '-', 'Color', [0, 0, 0]/255, 'LineWidth', 2, 'DisplayName', 'Mod. Advection fog'); hold on;

xlabel('\lambda (\mum)');
ylabel('Lidar Ratio (sr)');

xlim([wavelength(1), wavelength(end)]);
ylim([0, 50]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.01]);

export_fig(gcf, fullfile('D:\Coding\Matlab\LISMO\image\ensemble_sea_fog_optical_properties.png'), '-r300');