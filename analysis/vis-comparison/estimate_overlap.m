clc;
close all;

%% Parameter Definition
lidarFolder = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\lidar';
distOffset = 12 * 7.5;
filename = '20241101001000.pb';
olFile = 'overlap_20241101.mat';
linFitRange = [2000, 2500];
hRangeDisplay = [0, 3000];
savePath = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\vis-comparison';

% read lidar data
thisDate = datenum(filename(1:(end - 3)), 'yyyymmddHHMMSS');
data = readVisLidarL0(fullfile(lidarFolder, datestr(thisDate, 'yyyy-mm-dd'), filename));

%% Overlap Correction
range = transpose(((1:length(data.sig)) - 0.5) * data.hRes - distOffset);
bg = nanmean(data.sig((end - 300):end));
signal = data.sig - bg;
rcs = signal .* range.^2;
snr = sqrt(signal * 3000 * 50 * 1e-3) ./ sqrt(data.sig * 3000 * 50 * 1e-3);

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
ylim([1e3, 1e8]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
legend([p1, p2], 'Location', 'SouthEast');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'FontSize', 12);

subplot(212);
p1 = plot(range, smooth(rcs, 4) ./ smooth(fullRCS, 4), '-k', 'LineWidth', 2);
hold on;
p2 = plot([0, 1e5], [1, 1], '--r');

xlabel('距离 (米)');
ylabel('重叠因子');

xlim(hRangeDisplay);
ylim([-0.1, 1.2]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'FontSize', 12);

%% Save Overlap
ov = smooth(rcs, 4) ./ smooth(fullRCS, 4);
ov(range >= 1200) = 1;
save(olFile, 'ov', 'range');

if exist(savePath, 'dir')
    saveFile = fullfile(savePath, sprintf('lk_%s_overlap_factor.txt', datestr(data.startTime, 'yyyymmdd')));
    fid = fopen(saveFile, 'w');

    fprintf(fid, 'range (m) overlap_factor (smoothed by 4 range bins)\n');
    for iLine = 1:length(range)
        fprintf(fid, '%f %f\n', range(iLine), ov(iLine));
    end

    fclose(fid);
end