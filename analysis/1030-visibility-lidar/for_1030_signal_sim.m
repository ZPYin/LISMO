%------------------------------------------------------------------------------%
% Test the simulation of 1030 visibility lidar.
% 2024-03-24
%------------------------------------------------------------------------------%
close all;
global LISMO_VARS;

%% Parameter Definition
distArr = 7.5:15:20000;
eleAngle = 0;
laserWL = 1030;
pulseEn = 0.15;
eta = 0.07;
FWHMs = 3;
FOV = 0.4;
PB = 0.0;
darkCount = 300;
accShots = 5 * 10000;
visArr = 10000;
acqTime = 100;
optEffiEmit = -log10(0.8);
optEffiRecv = -log10(0.3);
saveFile = fullfile('C:\Users\ZPYin\Documents\INAST-Lidar-Group\能见度激光雷达\image\sig_1030_nighttime.png');

%% Signal Simulation
height = distArr * sin(eleAngle / 180 * pi);
tExt550 = vis2ext(visArr, 'method', 'mor') * ones(size(height));
mExt550 = 1.1288e-5 * ones(size(height));
mExt = mExt550 * (laserWL / 550) ^ -4;
aExt = (tExt550 - mExt550) * (laserWL / 550) ^ -0.8;
mBsc = mExt / (8 / 3 * pi);
aBsc = aExt / 64;
dataSim = LISMO_Model(distArr, 'tBsc', aBsc + mBsc, ...
                               'tExt', mExt + aExt, ...
                               'eleAngle', eleAngle, ...
                               'laserWL', laserWL, ...
                               'pulseEn', pulseEn, ...
                               'accShots', accShots, ...
                               'acqTime', acqTime, ...
                               'PB', PB, ...
                               'etaFR', eta, ...
                               'FOV_FR', FOV, ...
                               'NDEmit_FR', optEffiEmit, ...
                               'NDRecv_FR', optEffiRecv, ...
                               'FWHM', FWHMs, ...
                               'darkCountFR', darkCount, ...
                               'visible', 'off', ...
                               'ylim', [0, 4]);

%% Display
figure('Position', [0, 10, 400, 400], 'Units', 'Pixels', 'Color', 'w');

subfigs = subfigPos([0.15, 0.13, 0.82, 0.80], 2, 1, 0, 0.07);

% Signal
subplot('Position', subfigs(1, :), 'Units', 'Normalized');

hold on;
plot(dataSim.distArr / 1e3, (dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR) / (accShots * acqTime * 1e-3), 'LineWidth', 2);
plot([10, 10], [1e-10, 1e10], '-.r');
hold off;

xlim([0, 20]);
ylim([1e-5, 1e7]);

xlabel('');
ylabel('激光雷达信号 (MHz)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTick', 10.^(-5:3:7), 'YScale', 'log', 'Box', 'on', 'FontSize', 11);

% SNR
subplot('Position', subfigs(2, :), 'Units', 'Normalized');

hold on;
plot(dataSim.distArr / 1e3, (dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR) ./ sqrt(dataSim.N_FR_Poiss), 'LineWidth', 2);
plot(dataSim.distArr / 1e3, 3 * ones(size(dataSim.distArr)), '--k');
plot([10, 10], [1e-10, 1e10], '-.r');
hold off;

xlim([0, 20]);
ylim([1e-1, 1e5]);

xlabel('距离 (千米)');
ylabel('信噪比');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTick', 10.^(-1:2:5), 'YScale', 'log', 'Box', 'on', 'FontSize', 11);

export_fig(gcf, saveFile, '-r300');