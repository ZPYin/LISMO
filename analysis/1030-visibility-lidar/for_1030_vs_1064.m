%------------------------------------------------------------------------------%
% Compare the performance of 1030 and 1064 visibility lidar
% 2024-03-24
%------------------------------------------------------------------------------%
close all;
global LISMO_VARS;

%% Parameter Definition
distArr = 7.5:15:20000;
eleAngle = 0;
laserWL = [1030, 1064];
pulseEn = 0.15;
eta = [0.03, 0.07];
FWHMs = 3;
FOV = 0.4;
PB = 0.05;
darkCount = 300;
accShots = 5 * 10000;
visArr = 10000;
acqTime = 100;
optEffiEmit = -log10(0.8);
optEffiRecv = -log10(0.3);
saveFile = fullfile('C:\Users\ZPYin\Documents\INAST-Lidar-Group\�ܼ��ȼ����״�\image\sig_1030_vs1064_daytime.png');

%% Signal Simulation
dataSimAll = cell(0);
for iWL = 1:length(laserWL)
    height = distArr * sin(eleAngle / 180 * pi);
    tExt550 = vis2ext(visArr, 'method', 'mor') * ones(size(height));
    mExt550 = 1.1288e-5 * ones(size(height));
    mExt = mExt550 * (laserWL(iWL) / 550) ^ -4;
    aExt = (tExt550 - mExt550) * (laserWL(iWL) / 550) ^ -0.8;
    mBsc = mExt / (8 / 3 * pi);
    aBsc = aExt / 64;
    dataSim = LISMO_Model(distArr, 'tBsc', aBsc + mBsc, ...
                                'tExt', mExt + aExt, ...
                                'eleAngle', eleAngle, ...
                                'laserWL', laserWL(iWL), ...
                                'pulseEn', pulseEn, ...
                                'accShots', accShots, ...
                                'acqTime', acqTime, ...
                                'PB', PB, ...
                                'etaFR', eta(iWL), ...
                                'FOV_FR', FOV, ...
                                'NDEmit_FR', optEffiEmit, ...
                                'NDRecv_FR', optEffiRecv, ...
                                'FWHM', FWHMs, ...
                                'darkCountFR', darkCount, ...
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
p1 = plot(dataSimAll{1}.distArr / 1e3, (dataSimAll{1}.N_FR_Poiss - dataSimAll{1}.Nb_FR - dataSimAll{1}.Nd_FR) / (accShots * acqTime * 1e-3), 'color', 'm', 'LineWidth', 2, 'DisplayName', '1030');
p2 = plot(dataSimAll{2}.distArr / 1e3, (dataSimAll{2}.N_FR_Poiss - dataSimAll{2}.Nb_FR - dataSimAll{2}.Nd_FR) / (accShots * acqTime * 1e-3), 'color', 'b', 'LineWidth', 2, 'DisplayName', '1064');
plot([10, 10], [1e-10, 1e10], '-.r');
hold off;

xlim([0, 20]);
ylim([1e-5, 1e7]);

xlabel('');
ylabel('�����״��ź� (MHz)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTick', 10.^(-5:3:7), 'YScale', 'log', 'Box', 'on', 'FontSize', 11);
legend([p1, p2], 'Location', 'NorthEast');

% SNR
subplot('Position', subfigs(2, :), 'Units', 'Normalized');

hold on;
plot(dataSimAll{1}.distArr / 1e3, (dataSimAll{1}.N_FR_Poiss - dataSimAll{1}.Nb_FR - dataSimAll{1}.Nd_FR) ./ sqrt(dataSimAll{1}.N_FR_Poiss), 'color', 'm', 'LineWidth', 2, 'DisplayName', '1030');
plot(dataSimAll{2}.distArr / 1e3, (dataSimAll{2}.N_FR_Poiss - dataSimAll{2}.Nb_FR - dataSimAll{2}.Nd_FR) ./ sqrt(dataSimAll{2}.N_FR_Poiss), 'color', 'B', 'LineWidth', 2, 'DisplayName', '1064');
plot(dataSim.distArr / 1e3, 3 * ones(size(dataSim.distArr)), '--k');
plot([10, 10], [1e-10, 1e10], '-.r');
hold off;

xlim([0, 20]);
ylim([1e-1, 1e5]);

xlabel('���� (ǧ��)');
ylabel('�����');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTick', 10.^(-1:2:5), 'YScale', 'log', 'Box', 'on', 'FontSize', 11);

export_fig(gcf, saveFile, '-r300');