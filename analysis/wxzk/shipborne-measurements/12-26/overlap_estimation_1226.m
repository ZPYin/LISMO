clc;
close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\海试数据、效果图\12月26无人艇出海数据';
distOffset = 55.75;
tRange = [datenum(2023, 12, 26, 15, 0, 0), datenum(2023, 12, 26, 15, 2, 0)];   % temporal range for data input
olFile = 'overlap_20231226.mat';
linFitRange = [1300, 1800];
hRangeDisplay = [0, 3000];
savePath = '';

%% List Data Files
dataFilesSearched = listfile(dataFolder, '\w*.VIS', 1);
dataFiles = cell(0);
for iFile = 1:length(dataFilesSearched)
    fileBasename = basename(dataFilesSearched{iFile});
    fileTime = datenum(fileBasename(3:(end - 4)), 'yyyymmddHHMMSS');

    isInTRange = (fileTime >= tRange(1)) & (fileTime < tRange(2));
    if isInTRange
        dataFiles = cat(2, dataFiles, dataFilesSearched{iFile});
    end
end
data = readVIS(dataFiles);

%% Overlap Correction
range = ((1:data.nBins(1)) - 0.5) * data.hRes(1) - distOffset;
bg = nanmean(squeeze(data.rawSignal(:, 3, 2900:2980)), 2);
signal = squeeze(data.rawSignal(:, 3, :)) - repmat(bg, 1, data.nBins(1));
rcs = nansum(signal .* repmat(range, length(data.hRes), 1).^2, 1);
snr = (signal) ./ sqrt(squeeze(data.rawSignal(:, 3, :)));

%% Calculate RCS Slope
isInFit = (range >= linFitRange(1)) & (range <= linFitRange(2)) & (rcs > 0);
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

xlim(hRangeDisplay);
ylim([5e8, 1e10]);

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
    saveFile = fullfile(savePath, sprintf('wxzk_%s_overlap_factor.txt', datestr(data.mTime(1), 'yyyymmdd')));
    fid = fopen(saveFile, 'w');
    
    fprintf(fid, 'range (m) overlap_factor (smoothed by 4 range bins)\n');
    for iLine = 1:length(range)
        fprintf(fid, '%f %f\n', range(iLine), ov(iLine));
    end
    
    fclose(fid);
end