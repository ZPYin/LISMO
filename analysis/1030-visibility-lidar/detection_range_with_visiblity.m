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
pulseEn = 0.1;
eta = 0.07;
FWHMs = 3;
FOV = 0.4;
PB = [0.0, 0.05];
darkCount = 300;
accShots = 5 * 8000;
visArr = [500, 1000, 2000, 3000, 4000, 5000, 10000, 15000, 20000, 25000, 30000];
acqTime = 100;
optEffiEmit = -log10(0.8);
optEffiRecv = -log10(0.7);
saveFile = fullfile(LISMO_VARS.projectDir, 'image', 'sig_1030_detection_range_with_vis_night.png');

%% Signal Simulation
detectRange = NaN(length(PB), length(visArr));

for iPB = 1:length(PB)
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
                                    'PB', PB(iPB), ...
                                    'etaFR', eta, ...
                                    'FOV_FR', FOV, ...
                                    'NDEmit_FR', optEffiEmit, ...
                                    'NDRecv_FR', optEffiRecv, ...
                                    'FWHM', FWHMs, ...
                                    'darkCountFR', darkCount, ...
                                    'dTel', 0.075, ...
                                    'visible', 'off', ...
                                    'ylim', [0, 4]);

        snr = (dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR) ./ sqrt(dataSim.N_FR_Poiss);
        idx = find(snr <= 3, 1, 'first');
        detectRange(iPB, iVis) = dataSim.distArr(idx);
    end
end

%% detection range with visibility

% plot
figure('Position', [0, 10, 420, 230], 'Units', 'Pixels', 'Color', 'w');
hold on;
p1 = plot(visArr * 1e-3, detectRange(1, :) * 1e-3, 'color', 'b', 'LineWidth', 2, 'Marker', 's', 'MarkerFaceColor', 'b', 'DisplayName', '夜间');
p2 = plot(visArr * 1e-3, detectRange(2, :) * 1e-3, 'color', 'm', 'LineWidth', 2, 'Marker', 's', 'MarkerFaceColor', 'm', 'DisplayName', '白天');
plot([-1, 1e10], [10, 10], '--k');
hold off;

xlabel('能见度 (千米)');
ylabel('探测距离 (千米)');

xlim([0, 30]);
ylim([0, 20]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);
legend([p1, p2], 'Location', 'NorthWest');

export_fig(gcf, saveFile, '-r300');