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
marineSLCD = slwcDSD_mgd('N0', 74e6, 'gamma0', 8.6, 'Rm', 1.35, 'type', 'input');
continentalSLCD = slwcDSD_mgd('N0', 288e6, 'gamma0', 8.7, 'Rm', 0.65, 'type', 'input');

%% Cloud Optical Properties
marineSLBackscatter = NaN(length(wavelength), 1);
marineSLExtinction = NaN(length(wavelength), 1);
marineSLAsymFactor = NaN(length(wavelength), 1);
continentalSLBackscatter = NaN(length(wavelength), 1);
continentalSLExtinction = NaN(length(wavelength), 1);
continentalSLAsymFactor = NaN(length(wavelength), 1);

for iWL = 1:length(wavelength)
    marineSLBackscatter(iWL) = trapz(r, marineSLCD(r) .* qb(iWL, :) .* pi .* r.^2 * 1e-9 ) / (4 * pi);
    marineSLExtinction(iWL) = trapz(r, marineSLCD(r) .* qe(iWL, :) .* pi .* r.^2 * 1e-9);
    marineSLAsymFactor(iWL) = trapz(r, marineSLCD(r) .* qs(iWL, :) .* g(iWL, :) .* r.^2) ./ trapz(r, marineSLCD(r) .* qs(iWL, :) .* r.^2);

    continentalSLBackscatter(iWL) = trapz(r, continentalSLCD(r) .* qb(iWL, :) .* pi .* r.^2 * 1e-9) / (4 * pi);
    continentalSLExtinction(iWL) = trapz(r, continentalSLCD(r) .* qe(iWL, :) .* pi .* r.^2 * 1e-9);
    continentalSLAsymFactor(iWL) = trapz(r, continentalSLCD(r) .* qs(iWL, :) .* g(iWL, :) .* r.^2) ./ trapz(r, continentalSLCD(r) .* qs(iWL, :) .* r.^2);
end

%% liquid water content
lwc = trapz(r, marineSLCD(r) * 4/3 * pi .* r.^3) * 1e-15;   % g/m3
iwc = 5e-2;   % g/m3

% Display

[~, idx] = min(abs(wavelength - 0.532));
figure('Position', [0, 100, 500, 400], 'Units', 'Pixels', 'color', 'w');

hold on;
p1 = plot(10.^(-7:0.001:-2), lwc / (marineSLBackscatter(idx) .* 1e-3) .* 10.^(-7:0.001:-2), 'r-', 'LineWidth', 2, 'DisplayName', '液水云');
p2 = plot(10.^(-7:0.001:-2), iwc / (1e-2 / 25) .* 10.^(-7:0.001:-2), 'b-', 'LineWidth', 2, 'DisplayName', '冰云');
hold off;

xlim([1e-7, 1e-2]);
ylim([1e-8, 1]);

xlabel('(532nm) 后向散射系数 (m^{-1}sr^{-1})');
ylabel('冰/液水含量 (g/m^3)')

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'XTick', 10.^(-7:1:-2), 'YTick', 10.^(-8:2:0), 'Box', 'on', 'layer', 'top', 'TickLen', [0.02, 0.01], 'FontSize', 12, 'XScale', 'log', 'YScale', 'log');

text(0.05, 0.8, sprintf('最大穿透深度(COD=4)\n液水路径: %5.3f g m^{-2}\n冰水路径: %5.3f g m^{-2}', 1e-3 / (marineSLBackscatter(idx) .* 1e-3) * lwc * 160, 1e-3 / (1e-2 / 25) * iwc * 160), 'Units', 'normalized', 'FontSize', 12);
l = legend([p1, p2], 'Location', 'SouthEast');
l.FontSize = 12;