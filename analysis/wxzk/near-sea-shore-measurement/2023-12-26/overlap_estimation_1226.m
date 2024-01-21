clc;
close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\RawData\QingDao1\2023\12\26\H_SCAN15_135_2_20231226225006';
tRange = [datenum(2023, 12, 26,22, 50, 0), datenum(2023, 12, 26, 23, 0, 0)];   % temporal range for data input
olFile = 'overlap_20231226.mat';
distOffset = -48.75;
flagReadData = false;
linFitRange = [2000, 4000];
hRangeDisplay = [0, 5000];
savePath = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\wxzk\near-sea-shore-measurement\2023-12-26';
iCh = 1;

%% Read Data
if flagReadData
    lData = readVIS(dataFolder, 'tRange', tRange, 'isDir', true, 'debug', true);
end

%% Overlap Correction
range = ((1:lData.nBins(1)) - 0.5) * lData.hRes(1) + distOffset;
bg = nanmean(squeeze(lData.rawSignal(:, iCh, (end - 10):end)), 2);
signal = squeeze(lData.rawSignal(:, iCh, :)) - repmat(bg, 1, lData.nBins(1));
rcs = signal .* repmat(range, length(lData.hRes), 1).^2;
snr = (signal) ./ sqrt(squeeze(lData.rawSignal(:, iCh, :)));

%% Calculate RCS Slope
meanRCS = nanmean(rcs, 1);
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
ylim([1e8, 1e11]);

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
figure('Position', [0, 0, 600, 400], 'color', 'w', 'visible', 'on');

subplot('Position', [0.1, 0, 0.86, 0.86], 'Units', 'normalized');

rcs(snr <= 1) = NaN;
[~, p1] = polarPcolor(range / 1e3, lData.azimuthAng, rcs, 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 10], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 8e9], 'tickSize', 12, 'tickColor', 'm');
ylabel(p1, 'range cor. sig.');
set(p1, 'location', 'westoutside', 'FontSize', 12);
colormap(gca, myColormap('jetImage'));
text(0.3, 1.2, sprintf('%s', datestr(mean(lData.startTime), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'Bold');
text(0.6, 0, 'distance (km)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'Bold');

%% Save Overlap
ov = smooth(meanRCS, 4) ./ smooth(fullRCS, 4);
ov(range >= 2000) = 1;
save(olFile, 'ov', 'range');

if ~isempty(savePath)
    saveFile = fullfile(savePath, sprintf('wxzk_%s_overlap_factor.txt', datestr(lData.startTime(1), 'yyyymmdd')));
    fid = fopen(saveFile, 'w');
    
    fprintf(fid, 'range (m) overlap_factor (smoothed by 4 range bins)\n');
    for iLine = 1:length(range)
        fprintf(fid, '%f %f\n', range(iLine), ov(iLine));
    end
    
    fclose(fid);
end