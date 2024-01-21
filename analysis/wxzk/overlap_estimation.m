clc;
close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\RawData\QingDao1\2023\12\25\H_SCAN15_135_2_20231225210551';
distOffset = 48.75;
olFile = 'overlap_20231225.mat';
linFitRange = [2000, 2700];
hRangeDisplay = [0, 3000];
savePath = 'C:\Users\ZPYin\Desktop';

%% List Data Files
dataFiles = listfile(dataFolder, '\w*.VIS', 1);
dataFiles = dataFiles(1:30);
data = readVIS(dataFiles);

%% Overlap Correction
range = ((1:data.nBins(1)) - 0.5) * data.hRes(1) - distOffset;
bg = nanmean(squeeze(data.rawSignal(:, 1, 2900:2980)), 2);
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
title('2023年12月25日 00:17');

xlim(hRangeDisplay);
ylim([5e8, 1e12]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
legend([p1, p2], 'Location', 'NorthEast');

subplot(212);
p1 = plot(range, smooth(rcs, 4) ./ smooth(fullRCS, 4), '-k', 'LineWidth', 2);
hold on;
p2 = plot([0, 1e5], [1, 1], '--r');

xlabel('距离 (米)');
ylabel('重叠因子');

xlim(hRangeDisplay);
ylim([-0.1, 1.2]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

%% Save Overlap
ov = smooth(rcs, 4) ./ smooth(fullRCS, 4);
ov(range >= 2000) = 1;
save(olFile, 'ov', 'range');

if exist(savePath, 'dir')
    saveFile = fullfile(savePath, 'wxzk_20231225_overlap_factor.txt');
    fid = fopen(saveFile, 'w');
    
    fprintf(fid, 'range (m) overlap_factor (smoothed by 4 range bins)\n');
    for iLine = 1:length(range)
        fprintf(fid, '%f %f\n', range(iLine), ov(iLine));
    end
    
    fclose(fid);
end