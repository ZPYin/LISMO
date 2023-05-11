clc;
close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\RawData\QingDao\2023\4\17\H_SCAN0_120_2_20230417005424';
distOffset = 48.75;
linFitRange = [700, 1500];

%% List Data Files
dataFiles = listfile(dataFolder, '\w*.VIS', 1);
data = readVIS(dataFiles);

%% Overlap Correction
range = ((1:data.nBins(1)) - 0.5) * data.hRes(1) - distOffset;
bg = nanmean(squeeze(data.rawSignal(:, 1, 2900:2950)), 2);
signal = squeeze(data.rawSignal(:, 1, :)) - repmat(bg, 1, data.nBins(1));
rcs = nansum(signal .* repmat(range, length(data.hRes), 1).^2, 1);
snr = (signal) ./ sqrt(squeeze(data.rawSignal(:, 1, :)));

%% Calculate RCS Slope
isInFit = (range >= linFitRange(1)) & (range <= linFitRange(2));
fitRange = range(isInFit);
fitLogRCS = log(rcs(isInFit));
[offset, slope] = chi2fit(fitRange, fitLogRCS, log(sqrt(rcs(isInFit))));
positiveRCS = rcs;
positiveRCS(positiveRCS <= 0) = NaN;
fullRCS = exp(offset + slope * range);

%% Display
figure('Position', [0, 0, 700, 600], 'Units', 'Pixels', 'Color', 'w');

subplot(211);
p1 = semilogy(range, rcs, 'LineWidth', 2, 'DisplayName', '原始测量');
hold on;
p2 = semilogy(range, fullRCS, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', '重建后');

xlabel('距离 (米)');
ylabel('距离修正信号');
title('2023年4约17日 00:54扫描周期');

xlim([0, 2000]);
ylim([1e10, 1e12]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
legend([p1, p2], 'Location', 'NorthEast');

subplot(212);
p1 = plot(range, rcs ./ fullRCS, '-k', 'LineWidth', 2);
hold on;
p2 = plot([0, 1e5], [1, 1], '--r');

xlabel('距离 (米)');
ylabel('重叠因子');

xlim([0, 2000]);
ylim([-0.1, 1.2]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

%% Save Overlap
ov = rcs ./ fullRCS;
ov(range >= 600) = 1;
save('overlap.mat', 'ov', 'range');