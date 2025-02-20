% 对比中科光博反演的消光系数和能见度产品
% 作者：殷振平
% 日期：2025年2月9日

clc;
close all;

%% Parameter Definition
L1Path = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\ZKGB\results\L1';
zkgbDataPath = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\ZKGB\data\2025-01-01';
resPath = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\ZKGB\results\cmp';
figPath = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\ZKGB\results\cmp';
flagReadData = true;
distOffset = -15;

%% Single Scanning Comparison
scanPath = '20250101000758';

%% Read Data

% ZKGB extinction
zkgbExtFile = fullfile(zkgbDataPath, scanPath, sprintf('%s_1064Extinction.xlsx', scanPath));
thisMat = readmatrix(zkgbExtFile, 'range', [1 1]);
range = thisMat(1, 2:end) + distOffset;
thisZKGBExt = thisMat(2:end, 2:end);   % (km-1)

% ZKGB visibility
zkgbVisFile = fullfile(zkgbDataPath, scanPath, sprintf('%s_Visibility.xlsx', scanPath));
thisMat = readmatrix(zkgbVisFile, 'range', [1 1]);
thisZKGBVis = thisMat(2:end, 2:end);   % (km)

% ZKGB our retrieval
L1File = fullfile(L1Path, sprintf('%s_L1.mat', scanPath));
zkgbRet0 = load(L1File);

%% comparison
thisExtDev = (thisZKGBExt - zkgbRet0.extMat_Fernald * 1e3);
thisVisDev = (thisZKGBVis - zkgbRet0.visMat_Fernald * 1e-3);
thisExtRelDev = abs(thisExtDev ./ (zkgbRet0.extMat_Fernald * 1e3) * 100);
thisVisRelDev = abs(thisVisDev ./ (zkgbRet0.visMat_Fernald * 1e-3) * 100);

%% Display

% extinction
figure('Name', 'Cmp Extinction', 'Position', [100, 100, 1200, 400], 'color', 'w', 'visible', 'on');

subplot('Position', [0.02, 0.1, 0.23, 0.84], 'Units', 'normalized');

[~, p1] = polarPcolor(zkgbRet0.range / 1e3, transpose(zkgbRet0.angle), transpose(zkgbRet0.extMat_Fernald) * 1e3, 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 5], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 0.5], 'tickSize', 12, 'tickColor', 'm');
ylabel(p1, '(CMA) extinction (km-1)', 'FontSize', 12);
colormap(gca, myColormap('jetImage'));
set(p1, 'location', 'westoutside', 'FontSize', 12);

text(-0.2, -0.09, sprintf('%s', datestr(zkgbRet0.mTime(1), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'Bold');

subplot('Position', [0.27, 0.1, 0.23, 0.84], 'Units', 'normalized');

[~, p1] = polarPcolor(zkgbRet0.range / 1e3, transpose(zkgbRet0.angle), transpose(thisZKGBExt), 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 5], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 0.5], 'tickSize', 12, 'tickColor', 'm');
ylabel(p1, '(ZKGB) extinction (km-1)', 'FontSize', 12);
colormap(gca, myColormap('jetImage'));
set(p1, 'location', 'westoutside', 'FontSize', 12);

subplot('Position', [0.52, 0.1, 0.23, 0.84], 'Units', 'normalized');

[~, p1] = polarPcolor(zkgbRet0.range / 1e3, transpose(zkgbRet0.angle), transpose(thisExtDev), 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 5], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [-0.1, 0.1], 'tickSize', 12, 'tickColor', 'm');
ylabel(p1, '\Delta extinction (km-1)', 'FontSize', 12);
colormap(gca, myColormap('jetImage'));
set(p1, 'location', 'westoutside', 'FontSize', 12);

subplot('Position', [0.77, 0.1, 0.23, 0.84], 'Units', 'normalized');

[~, p1] = polarPcolor(zkgbRet0.range / 1e3, transpose(zkgbRet0.angle), transpose(thisExtRelDev), 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 5], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 30], 'tickSize', 12, 'tickColor', 'm');
ylabel(p1, 'Rel.Dev. (%%)', 'FontSize', 12);
colormap(gca, myColormap('jetImage'));
set(p1, 'location', 'westoutside', 'FontSize', 12);

export_fig(gcf, fullfile(figPath, sprintf('%s_cmp_Extinction.png', scanPath)), '-r300');

% visibility
figure('Name', 'Cmp Visibility', 'Position', [100, 100, 1200, 400], 'color', 'w', 'visible', 'on');

subplot('Position', [0.02, 0.1, 0.23, 0.84], 'Units', 'normalized');

[~, p1] = polarPcolor(zkgbRet0.range / 1e3, transpose(zkgbRet0.angle), transpose(zkgbRet0.visMat_Fernald) * 1e-3, 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 5], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 30], 'tickSize', 12, 'tickColor', 'm');
ylabel(p1, '(CMA) visibility (km)', 'FontSize', 12);
colormap(gca, myColormap('jetImage'));
set(p1, 'location', 'westoutside', 'FontSize', 12);

text(-0.2, -0.09, sprintf('%s', datestr(zkgbRet0.mTime(1), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'Bold');

subplot('Position', [0.27, 0.1, 0.23, 0.84], 'Units', 'normalized');

[~, p1] = polarPcolor(zkgbRet0.range / 1e3, transpose(zkgbRet0.angle), transpose(thisZKGBVis), 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 5], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 30], 'tickSize', 12, 'tickColor', 'm');
ylabel(p1, '(ZKGB) visibility (km)', 'FontSize', 12);
colormap(gca, myColormap('jetImage'));
set(p1, 'location', 'westoutside', 'FontSize', 12);

subplot('Position', [0.52, 0.1, 0.23, 0.84], 'Units', 'normalized');

[~, p1] = polarPcolor(zkgbRet0.range / 1e3, transpose(zkgbRet0.angle), transpose(thisVisDev), 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 5], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [-3, 3], 'tickSize', 12, 'tickColor', 'm');
ylabel(p1, '\Delta visibility (km)', 'FontSize', 12);
colormap(gca, myColormap('jetImage'));
set(p1, 'location', 'westoutside', 'FontSize', 12);

subplot('Position', [0.77, 0.1, 0.23, 0.84], 'Units', 'normalized');

[~, p1] = polarPcolor(zkgbRet0.range / 1e3, transpose(zkgbRet0.angle), transpose(thisVisRelDev), 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 5], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 30], 'tickSize', 12, 'tickColor', 'm');
ylabel(p1, 'Rel.Dev. (%%)', 'FontSize', 12);
colormap(gca, myColormap('jetImage'));
set(p1, 'location', 'westoutside', 'FontSize', 12);

export_fig(gcf, fullfile(figPath, sprintf('%s_cmp_Visibility.png', scanPath)), '-r300');

%% Full Day Comparison
iPrf = 1;

%% Iteration
if flagReadData
    scanPaths = listdir(zkgbDataPath, '\w{14}', 1);
    nBins = 2500;
    mTime = NaN(1, length(scanPaths));
    fullDayZKGBExt0 = NaN(length(scanPaths), nBins);
    fullDayZKGBVis0 = NaN(length(scanPaths), nBins);
    fullDayZKGBExt1 = NaN(length(scanPaths), nBins);
    fullDayZKGBVis1 = NaN(length(scanPaths), nBins);
    for iScan = 1:length(scanPaths)
        fprintf('Finished %6.2f%%: %s...\n', iScan / (length(scanPaths) - 1) * 100, basename(scanPaths{iScan}));

        thisScanTimeStr = basename(scanPaths{iScan});
        zkgbExtFile1 = fullfile(scanPaths{iScan}, sprintf('%s_1064Extinction.xlsx', thisScanTimeStr));
        zkgbVisFile1 = fullfile(scanPaths{iScan}, sprintf('%s_Visibility.xlsx', thisScanTimeStr));
        zkgbFile0 = fullfile(L1Path, sprintf('%s_L1.mat', thisScanTimeStr));

        % ZKGB extinction
        thisMat = readmatrix(zkgbExtFile1, 'range', [1 1]);
        thisZKGBExt1 = thisMat(2:end, 2:end);   % (km-1)

        % ZKGB visibility
        thisMat = readmatrix(zkgbVisFile1, 'range', [1 1]);
        thisZKGBVis1 = thisMat(2:end, 2:end);   % (km)

        % ZKGB our retrieval
        zkgbRet0 = load(zkgbFile0);

        range = zkgbRet0.range(1:nBins);
        mTime(iScan) = zkgbRet0.mTime(iPrf);
        fullDayZKGBExt0(iScan, :) = zkgbRet0.extMat_Fernald(iPrf, 1:nBins) * 1e3;
        fullDayZKGBVis0(iScan, :) = zkgbRet0.visMat_Fernald(iPrf, 1:nBins) * 1e-3;
        fullDayZKGBExt1(iScan, :) = thisZKGBExt1(iPrf, 1:nBins);
        fullDayZKGBVis1(iScan, :) = thisZKGBVis1(iPrf, 1:nBins);
    end
end

%% Comparison
fullDayExtDev = (fullDayZKGBExt1 - fullDayZKGBExt0);
fullDayVisDev = (fullDayZKGBVis1 - fullDayZKGBVis0);
fullDayExtRelDev = abs(fullDayExtDev ./ (fullDayZKGBExt0) * 100);
fullDayVisRelDev = abs(fullDayVisDev ./ (fullDayZKGBVis0) * 100);

% some stats.
devStatsRange = [500, 3000];
isInStatsRange = (range >= devStatsRange(1)) & (range <= devStatsRange(2));
meanExtDev = mean(fullDayExtDev(:, isInStatsRange), 2);
meanVisDev = mean(fullDayVisDev(:, isInStatsRange), 2);
meanExtRelDev = mean(fullDayExtRelDev(:, isInStatsRange), 2);
meanVisRelDev = mean(fullDayVisRelDev(:, isInStatsRange), 2);

%% Output
fprintf('Some statistics....\n');
fprintf('Mean Extinction Dev. %6.2f km-1 (%d-%d m)\n', mean(meanExtDev), devStatsRange(1), devStatsRange(2));
fprintf('Mean Visibility Dev. %6.2f km (%d-%d m)\n', mean(meanVisDev), devStatsRange(1), devStatsRange(2));
fprintf('Mean Extinction Rel. Dev. %6.2f%% (%d-%d m)\n', mean(meanExtRelDev), devStatsRange(1), devStatsRange(2));
fprintf('Mean Visibility Rel. Dev. %6.2f%% (%d-%d m)\n', mean(meanVisRelDev), devStatsRange(1), devStatsRange(2));
save(fullfile(resPath, 'fullDayCmp.mat'), 'range', 'mTime', 'fullDayZKGBExt0', 'fullDayZKGBVis0', 'fullDayZKGBExt1', 'fullDayZKGBVis1', 'fullDayExtDev', 'fullDayVisDev', 'fullDayExtRelDev', 'fullDayVisRelDev', 'iPrf');

%% Display

% Mean extinction and visiblity dev
figure('Name', 'Mean Ext. Dev', 'Position', [100, 100, 450, 300], 'Units', 'Pixels', 'Color', 'w', 'visible', 'on');

hold on;
plot(mTime, meanExtDev, '-b', 'LineWidth', 1);
plot([0, 1e10], [0, 0], '--k');
hold off;

xlabel('Time');
ylabel('Mean Ext. Dev. (km-1)');
title('Mean Ext. Dev.');

xlim([min(mTime(1)), max(mTime(end))]);
ylim([-0.1, 0.1]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'FontSize', 12, 'Box', 'on');

datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

export_fig(gcf, fullfile(figPath, 'meanExtDev.png'), '-r300');

figure('Name', 'Mean Vis. Dev', 'Position', [100, 100, 450, 300], 'Units', 'Pixels', 'Color', 'w', 'visible', 'on');

hold on;
plot(mTime, meanVisDev, '-b', 'LineWidth', 1);
plot([0, 1e10], [0, 0], '--k');
hold off;

xlabel('Time');
ylabel('Mean Vis. Dev. (km)');
title('Mean Vis. Dev.');

xlim([min(mTime(1)), max(mTime(end))]);
ylim([-10, 10]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'FontSize', 12, 'Box', 'on');

datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

export_fig(gcf, fullfile(figPath, 'meanVisDev.png'), '-r300');

figure('Name', 'Mean Ext. Rel. Dev.', 'Position', [100, 100, 450, 300], 'Units', 'Pixels', 'Color', 'w', 'visible', 'on');

hold on;
plot(mTime, meanExtRelDev, '-b', 'LineWidth', 1);
plot([0, 1e10], [0, 0], '--k');
hold off;

xlabel('Time');
ylabel('Mean Ext. Rel. Dev. (%)');
title('Mean Ext. Rel. Dev.');

xlim([min(mTime(1)), max(mTime(end))]);
ylim([0, 50]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'FontSize', 12, 'Box', 'on');

datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

export_fig(gcf, fullfile(figPath, 'meanExtRelDev.png'), '-r300');

figure('Name', 'Mean Vis. Rel. Dev.', 'Position', [100, 100, 450, 300], 'Units', 'Pixels', 'Color', 'w', 'visible', 'on');

hold on;
plot(mTime, meanVisRelDev, '-b', 'LineWidth', 1);
plot([0, 1e10], [0, 0], '--k');
hold off;

xlabel('Time');
ylabel('Mean Vis. Rel. Dev. (%)');
title('Mean Vis. Rel. Dev.');

xlim([min(mTime(1)), max(mTime(end))]);
ylim([0, 50]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'FontSize', 12, 'Box', 'on');

datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

export_fig(gcf, fullfile(figPath, 'meanVisRelDev.png'), '-r300');

% Extinction
figure('Name', 'Full Day Comparison Extinction', 'Position', [100, 100, 800, 400], 'Units', 'Pixels', 'color', 'w', 'visible', 'on');

subFigs = subfigPos([0.06, 0.1, 0.8, 0.83], 2, 2, 0.1, 0.07);

subplot('Position', subFigs(1, :), 'Units', 'Normalized');

hold on;
p1 = pcolor(mTime, range(1:nBins) * 1e-3, transpose(fullDayZKGBExt0));
hold off;
p1.EdgeColor = 'none';
caxis([0, 0.5]);

xlim([min(mTime(1)), max(mTime(end))]);
ylim([0, 5]);

colormap('jet');
% xlabel('Time');
ylabel('Range (km)');
title('(CMA) Extinction (km-1)');

set(gca, 'XMinorTick', 'on', 'XTickLabel', '', 'YMinorTick', 'on', 'FontSize', 12, 'TickDir', 'out', 'Box', 'on', 'layer', 'top');

% datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

colorbar('Position', [subFigs(1, 1) + subFigs(1, 3) + 0.01, subFigs(1, 2), 0.01, subFigs(1, 4)], 'Units', 'Normalized');

subplot('Position', subFigs(2, :), 'Units', 'Normalized');

hold on;
p1 = pcolor(mTime, range(1:nBins) * 1e-3, transpose(fullDayZKGBExt1));
hold off;
p1.EdgeColor = 'none';
caxis([0, 0.5]);

xlim([min(mTime(1)), max(mTime(end))]);
ylim([0, 5]);

colormap('jet');
% xlabel('Time');
ylabel('Range (km)');
title('(ZKGB) Extinction (km-1)');

set(gca, 'XMinorTick', 'on', 'XTickLabel', '', 'YMinorTick', 'on', 'FontSize', 12, 'TickDir', 'out', 'Box', 'on', 'layer', 'top');

% datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

colorbar('Position', [subFigs(2, 1) + subFigs(2, 3) + 0.01, subFigs(2, 2), 0.01, subFigs(2, 4)], 'Units', 'Normalized');

subplot('Position', subFigs(3, :), 'Units', 'Normalized');

hold on;
p1 = pcolor(mTime, range(1:nBins) * 1e-3, transpose(fullDayExtDev));
hold off;
p1.EdgeColor = 'none';
caxis([-0.1, 0.1]);

xlim([min(mTime(1)), max(mTime(end))]);
ylim([0, 5]);

colormap('jet');
% xlabel('Time');
ylabel('Range (km)');
title('\Delta Extinction (km-1)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'FontSize', 12, 'TickDir', 'out', 'Box', 'on', 'layer', 'top');

datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

colorbar('Position', [subFigs(3, 1) + subFigs(3, 3) + 0.01, subFigs(3, 2), 0.01, subFigs(3, 4)], 'Units', 'Normalized');

subplot('Position', subFigs(4, :), 'Units', 'Normalized');

hold on;
p1 = pcolor(mTime, range(1:nBins) * 1e-3, transpose(fullDayExtRelDev));
hold off;
p1.EdgeColor = 'none';
caxis([0, 30]);

xlim([min(mTime(1)), max(mTime(end))]);
ylim([0, 5]);

colormap('jet');
% xlabel('Time');
ylabel('Range (km)');
title('Rel.Dev. (%)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'FontSize', 12, 'TickDir', 'out', 'Box', 'on', 'layer', 'top');

datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

colorbar('Position', [subFigs(4, 1) + subFigs(4, 3) + 0.01, subFigs(4, 2), 0.01, subFigs(4, 4)], 'Units', 'Normalized');

export_fig(gcf, fullfile(figPath, 'fullDayCmp_extinction.png'), '-r300');

% Visibility
figure('Name', 'Full Day Comparison Visibility', 'Position', [100, 100, 800, 400], 'Units', 'Pixels', 'color', 'w', 'visible', 'on');

subFigs = subfigPos([0.06, 0.1, 0.8, 0.83], 2, 2, 0.1, 0.07);

subplot('Position', subFigs(1, :), 'Units', 'Normalized');

hold on;
p1 = pcolor(mTime, range(1:nBins) * 1e-3, transpose(fullDayZKGBVis0));
hold off;
p1.EdgeColor = 'none';
caxis([0, 30]);

xlim([min(mTime(1)), max(mTime(end))]);
ylim([0, 5]);

colormap('jet');
% xlabel('Time');
ylabel('Range (km)');
title('(CMA) Visibility (km)');

set(gca, 'XMinorTick', 'on', 'XTickLabel', '', 'YMinorTick', 'on', 'FontSize', 12, 'TickDir', 'out', 'Box', 'on', 'layer', 'top');

% datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

colorbar('Position', [subFigs(1, 1) + subFigs(1, 3) + 0.01, subFigs(1, 2), 0.01, subFigs(1, 4)], 'Units', 'Normalized');

subplot('Position', subFigs(2, :), 'Units', 'Normalized');

hold on;
p1 = pcolor(mTime, range(1:nBins) * 1e-3, transpose(fullDayZKGBVis1));
hold off;
p1.EdgeColor = 'none';
caxis([0, 30]);

xlim([min(mTime(1)), max(mTime(end))]);
ylim([0, 5]);

colormap('jet');
% xlabel('Time');
ylabel('Range (km)');
title('(ZKGB) Visibility (km)');

set(gca, 'XMinorTick', 'on', 'XTickLabel', '', 'YMinorTick', 'on', 'FontSize', 12, 'TickDir', 'out', 'Box', 'on', 'layer', 'top');

% datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

colorbar('Position', [subFigs(2, 1) + subFigs(2, 3) + 0.01, subFigs(2, 2), 0.01, subFigs(2, 4)], 'Units', 'Normalized');

subplot('Position', subFigs(3, :), 'Units', 'Normalized');

hold on;
p1 = pcolor(mTime, range(1:nBins) * 1e-3, transpose(fullDayVisDev));
hold off;
p1.EdgeColor = 'none';
caxis([-3. 3]);

xlim([min(mTime(1)), max(mTime(end))]);
ylim([0, 5]);

colormap('jet');
% xlabel('Time');
ylabel('Range (km)');
title('\Delta Visibility (km)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'FontSize', 12, 'TickDir', 'out', 'Box', 'on', 'layer', 'top');

datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

colorbar('Position', [subFigs(3, 1) + subFigs(3, 3) + 0.01, subFigs(3, 2), 0.01, subFigs(3, 4)], 'Units', 'Normalized');

subplot('Position', subFigs(4, :), 'Units', 'Normalized');

hold on;
p1 = pcolor(mTime, range(1:nBins) * 1e-3, transpose(fullDayVisRelDev));
hold off;
p1.EdgeColor = 'none';
caxis([0, 30]);

xlim([min(mTime(1)), max(mTime(end))]);
ylim([0, 5]);

colormap('jet');
% xlabel('Time');
ylabel('Range (km)');
title('Rel.Dev. (%)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'FontSize', 12, 'TickDir', 'out', 'Box', 'on', 'layer', 'top');

datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

colorbar('Position', [subFigs(4, 1) + subFigs(4, 3) + 0.01, subFigs(4, 2), 0.01, subFigs(4, 4)], 'Units', 'Normalized');

export_fig(gcf, fullfile(figPath, 'fullDayCmp_visibility.png'), '-r300');