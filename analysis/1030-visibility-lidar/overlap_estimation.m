clc;
close all;

%% Parameter Definition
dataFile = 'G:\backup\vis-lidar\20260428-tianjin-evaluation\前散仪对比数据\VIS2\L0\2026-05-02_vis_lidar_l0.mat';
olFile = 'overlap_20260425.mat';
linFitRange = [800, 1000];
prfIdx = 720:730;
hRangeDisplay = [0, 2000];
savePath = 'C:\Users\zhenp\Documents\Coding\Matlab\LISMO\analysis\1030-visibility-lidar';
olFile = 'VIS2_overlap_factor.txt';

%% List Data Files
data = load(dataFile);

%% Preprocessing
range = data.range;
bg = nanmean(data.lidarSig(1500:1550, :), 1);
signal = data.lidarSig - repmat(bg, length(data.range), 1);
rcsAll = signal .* repmat(range, 1, length(data.mTime)).^2;
rcs = sum(rcsAll(:, prfIdx), 2);

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
title(sprintf('%s', datestr(mean(data.mTime(prfIdx)), 'yyyy-mm-dd HH:MM')));

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
    saveFile = fullfile(savePath, olFile);
    fid = fopen(saveFile, 'w');
    
    fprintf(fid, 'range (m) overlap_factor (smoothed by 4 range bins)\n');
    for iLine = 1:length(range)
        fprintf(fid, '%f %f\n', range(iLine), ov(iLine));
    end
    
    fclose(fid);
end