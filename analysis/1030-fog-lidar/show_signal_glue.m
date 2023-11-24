clc;
close all;
global LISMO_VARS;

%% 参数设置
dataFolder = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\data\1030\2023-11-07';
deadtime = 22;
firstRangeBin = 15;
tRange = [datenum(2023, 11, 07, 12, 0, 0), datenum(2023, 11, 07, 12, 8, 0)];
linearFitRange = [1200, 3000];

%% 读取数据
data = readALADats(dataFolder, 'tRange', tRange, 'nMaxBin', 2100);

%% 数据预处理
rawSigPCR = data.rawSignal / (50 * data.nShots(1) * 1e-3);
corSigPCR = rawSigPCR ./ (1 - deadtime * rawSigPCR * 1e-3);
corSigPC = corSigPCR * 50 * data.nShots(1) * 1e-3;
bg = nanmean(corSigPC(:, (end - 80):(end - 10), :), 2);
sigPC = corSigPC - repmat(bg, 1, size(corSigPC, 2), 1);
range = ((1:size(sigPC, 2)) - firstRangeBin) * data.hRes(1);
rcs = sigPC .* repmat(reshape(range, 1, length(range), 1), size(sigPC, 1), 1, size(sigPC, 3)).^2;

%% 线性拟合信号重构
meanRCS = mean(rcs, 3);
positiveRCS = squeeze(meanRCS(2, :));
positiveRCS(positiveRCS <= 0) = NaN;
isInFitRange = (range < linearFitRange(2)) & (range > linearFitRange(1));
fitRange = range(isInFitRange);
fitLogRCS = log(positiveRCS(isInFitRange));
[offset, slope] = chi2fit(fitRange, fitLogRCS, log(sqrt(positiveRCS(isInFitRange))));
fullRCS = exp(offset + slope * range);

%% 重叠因子计算
meanRCSFar = squeeze(nanmean(rcs(2, :, :), 3));
meanRCSNear = squeeze(nanmean(rcs(1, :, :), 3));
[maxRCSNear, maxRCSIdx] = max(meanRCSNear(1:100));
ratioFar = (nansum(fullRCS(isInFitRange)) / nansum(meanRCSFar(isInFitRange)));
overlapFar = smooth(meanRCSFar, 8) ./ smooth(fullRCS, 8) * ratioFar;
ratioNear = maxRCSNear / fullRCS(maxRCSIdx);
overlapNear = meanRCSNear ./ fullRCS * 1 / ratioNear * 0.95;

%% 显示
figure('Position', [0, 30, 500, 500], 'Units', 'Pixels', 'Color', 'w');

subfig = subfigPos([0.12, 0.1, 0.84, 0.8], 2, 1, 0, 0.1);

subplot('Position', subfig(1, :), 'Units', 'Pixels');
hold on;
p1 = plot(range, meanRCSFar / ratioFar, 'Color', 'r', 'LineWidth', 2, 'DisplayName', '远场信号');
p2 = plot(range, meanRCSNear / ratioNear, 'Color', 'b', 'LineWidth', 2, 'DisplayName', sprintf('近场信号x%d', 12));
p3 = plot(range, fullRCS, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', '理论信号');

xlabel('距离 (米)');
ylabel('距离修正信号 (a.u.)');
title('激光测雾雷达信号与重叠因子');

xlim([0, 2000]);
ylim([1e5, 1e12]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log');
legend([p1, p2, p3], 'Location', 'southeast');

subplot('Position', subfig(2, :), 'Units', 'Pixels');

hold on;
p1 = plot(range, overlapFar, 'color', 'r', 'LineWidth', 2, 'DisplayName', '重叠因子');
p2 = plot(range, overlapNear, 'color', 'b', 'LineWidth', 2, 'DisplayName', '重叠因子');
p3 = plot(range(18), overlapNear(18), 'color', 'b', 'Marker', '.', 'MarkerSize', 20);
plot([0, 1e10], [1, 1], '--k');
plot([0, 1e10], [0.1, 0.1], '--k');
plot([range(18), range(18)], [-1, 10], '--k');
text(range(18), overlapNear(18) + 0.3, sprintf('%d m', range(18)), 'Units', 'data');

xlabel('距离 (米)');
ylabel('重叠因子 (a.u.)');

xlim([0, 2000]);
ylim([-0.05, 1.2]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

export_fig(gcf, fullfile(sprintf('%s-overlap-function.pdf', datestr(data.mTime(1), 'yyyy-mm-dd'))));
export_fig(gcf, fullfile(sprintf('%s-overlap-function.png', datestr(data.mTime(1), 'yyyy-mm-dd'))), '-r300');
save(fullfile(LISMO_VARS.projectDir, 'data', sprintf('1030-fog-lidar-overlap-%s.mat', datestr(data.mTime(1), 'yyyy-mm-dd'))), 'range', 'overlapFar', 'overlapNear');