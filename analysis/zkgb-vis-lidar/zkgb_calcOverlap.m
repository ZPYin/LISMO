% 计算中科光博能见度激光雷达重叠因子
% 作者：殷振平
% 时间：2025-02-08

clc;
close all;

%% Parameter Definition
dataFile = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\ZKGB\data\2025-01-01\20250101192740\20250101192740_1064Signal.xlsx';
iPrf = 60:85;
olFile = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\zkgb-vis-lidar\overlap_20250101.mat';
distOffset = -15;
flagReadData = true;
linFitRange = [400, 700];
hRangeDisplay = [0, 5000];
savePath = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\ZKGB\results';

%% Read Data
if flagReadData
    lData = readVIS_ZKGB(dataFile, 'debug', false);
end

%% Overlap Correction
range = lData.range + distOffset;
bg = nanmean(lData.signal(:, (end - 20):end), 2);
noise = nanstd(lData.signal(:, (end - 20):end), 2);
signal = lData.signal - repmat(bg, 1, lData.nBins);
rcs = signal .* repmat(range, size(lData.signal, 1), 1).^2;
snr = signal ./ repmat(noise, 1, lData.nBins);

%% Calculate RCS Slope
meanRCS = nanmean(rcs(iPrf, :), 1);
isInFit = (range >= linFitRange(1)) & (range <= linFitRange(2)) & (meanRCS > 0);
fitRange = range(isInFit);
fitLogRCS = log(meanRCS(isInFit));
[offset, slope] = chi2fit(fitRange, fitLogRCS, log(sqrt(rcs(isInFit))));
positiveRCS = meanRCS;
positiveRCS(positiveRCS <= 0) = NaN;
fullRCS = exp(offset + slope * range);

%% Display

% overlap function
figure('Position', [0, 0, 600, 500], 'Units', 'Pixels', 'Color', 'w');

subplot(211);
p1 = semilogy(range, meanRCS, 'LineWidth', 2, 'DisplayName', '原始测量');
hold on;
p2 = semilogy(range, fullRCS, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', '重建后');

xlabel('距离 (米)');
ylabel('距离修正信号');

xlim(hRangeDisplay);
ylim([1e5, 1e7]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
legend([p1, p2], 'Location', 'NorthEast');

subplot(212);
plot(range, smooth(meanRCS, 4) ./ smooth(fullRCS, 4), '-k', 'LineWidth', 2);
hold on;
plot([0, 1e5], [1, 1], '--r');

xlabel('距离 (米)');
ylabel('重叠因子');

xlim(hRangeDisplay);
ylim([-0.1, 1.2]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

% measurements
figure('Position', [100, 100, 600, 400], 'color', 'w', 'visible', 'on');

subplot('Position', [0.1, 0.1, 0.86, 0.86], 'Units', 'normalized');

% rcs(snr <= 1) = NaN;
[~, p1] = polarPcolor(range / 1e3, transpose(lData.angle), transpose(rcs), 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 3], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 0.6e7], 'tickSize', 12, 'tickColor', 'm');
ylabel(p1, 'range cor. sig.');
set(p1, 'location', 'westoutside', 'FontSize', 12);
colormap(gca, myColormap('jetImage'));
text(0.3, 1.2, sprintf('%s', datestr(mean(lData.mTime), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'Bold');
text(0.6, -0.06, 'distance (km)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'Bold');

%% Save Overlap
ov = smooth(meanRCS, 4) ./ smooth(fullRCS, 4);
ov(range >= 500) = 1;
save(olFile, 'ov', 'range');

if ~isempty(savePath)
    saveFile = fullfile(savePath, sprintf('zkgb_%s_overlap_factor.txt', datestr(lData.mTime(1), 'yyyymmdd')));
    fid = fopen(saveFile, 'w');
    
    fprintf(fid, 'range (m) overlap_factor (smoothed by 4 range bins)\n');
    for iLine = 1:length(range)
        fprintf(fid, '%f %f\n', range(iLine), ov(iLine));
    end
    
    fclose(fid);
end