%------------------------------------------------------------------------------%
% Simulate the detection range under different visibilities.
% 2024-03-24
%------------------------------------------------------------------------------%
close all;
global LISMO_VARS;

%% Parameter Definition
distArr = 7.5:15:50000;
eleAngle = 0;
laserWL = 1030;
pulseEn = 0.15;
eta = 0.07;
FWHMs = 3;
FOV = 0.4;
PB = 0.05;
darkCount = 300;
accShots = 60 * 8000;
visArr = [500, 1000, 2000, 3000, 4000, 5000, 10000, 15000, 20000, 25000, 30000, 40000, 50000];
acqTime = 100;
optEffiEmit = -log10(0.8);
optEffiRecv = -log10(0.3);
saveFile = fullfile(LISMO_VARS.projectDir, 'image', 'sig_1030_detection_range_with_vis.png');

%% Signal Simulation
dataSimAll = cell(0);
for iVis = 1:length(visArr)
    height = distArr * sin(eleAngle / 180 * pi);
    tExt550 = vis2ext(visArr(iVis), 'method', 'mor') * ones(size(height));
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
                                'dTel', 0.175, ...
                                'visible', 'off', ...
                                'ylim', [0, 4]);
    dataSimAll = cat(2, dataSimAll, dataSim);
end

%% Display
figure('Position', [0, 10, 400, 400], 'Units', 'Pixels', 'Color', 'w');

subfigs = subfigPos([0.15, 0.13, 0.82, 0.80], 2, 1, 0, 0.07);

% Signal
subplot('Position', subfigs(1, :), 'Units', 'Normalized');

hold on;
p1 = plot(dataSimAll{1}.distArr / 1e3, (dataSimAll{1}.N_FR_Poiss) / (accShots * acqTime * 1e-3), 'Color', 'b', 'LineWidth', 2, 'DisplayName', '夜间');
p2 = plot(dataSimAll{2}.distArr / 1e3, (dataSimAll{2}.N_FR_Poiss) / (accShots * acqTime * 1e-3), 'Color', 'm', 'LineWidth', 2, 'DisplayName', '白天');
plot([10, 10], [1e-10, 1e10], '-.r');
hold off;

xlim([0, 20]);
ylim([1e-5, 1e7]);

xlabel('');
ylabel('激光雷达信号 (MHz)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTick', 10.^(-5:3:7), 'YScale', 'log', 'Box', 'on', 'FontSize', 11);

legend([p1, p2], 'Location', 'NorthEast');

% SNR
subplot('Position', subfigs(2, :), 'Units', 'Normalized');

hold on;
plot(dataSimAll{1}.distArr / 1e3, (dataSimAll{1}.N_FR_Poiss - dataSimAll{1}.Nb_FR - dataSimAll{1}.Nd_FR) ./ sqrt(dataSimAll{1}.N_FR_Poiss), 'Color', 'b', 'LineWidth', 2);
plot(dataSimAll{2}.distArr / 1e3, (dataSimAll{2}.N_FR_Poiss - dataSimAll{2}.Nb_FR - dataSimAll{2}.Nd_FR) ./ sqrt(dataSimAll{2}.N_FR_Poiss), 'Color', 'm', 'LineWidth', 2);
plot(dataSimAll{1}.distArr / 1e3, 3 * ones(size(dataSimAll{1}.distArr)), '--k');
plot([10, 10], [1e-10, 1e10], '-.r');
hold off;

xlim([0, 20]);
ylim([1e-1, 1e5]);

xlabel('距离 (千米)');
ylabel('信噪比');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTick', 10.^(-1:2:5), 'YScale', 'log', 'Box', 'on', 'FontSize', 11);

%% detection range with visibility
detectRange = zeros(size(visArr));
for iVis = 1:length(visArr)
    snr = (dataSimAll{iVis}.N_FR_Poiss - dataSimAll{iVis}.Nb_FR - dataSimAll{iVis}.Nd_FR) ./ sqrt(dataSimAll{iVis}.N_FR_Poiss);
    idx = find(snr <= 3, 1, 'first');
    detectRange(iVis) = dataSimAll{iVis}.distArr(idx);
end

% plot
figure('Position', [0, 10, 420, 230], 'Units', 'Pixels', 'Color', 'w');
hold on;
plot(visArr * 1e-3, detectRange * 1e-3, 'LineWidth', 2, 'Marker', 's', 'MarkerFaceColor', 'b');
hold off;

xlabel('能见度 (千米)');
ylabel('探测距离 (千米)');

xlim([0, 50]);
ylim([0, 15]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

export_fig(gcf, saveFile, '-r300');