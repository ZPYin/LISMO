% 分析给定时间范围内的岸基扫描能见度雷达数据，并将结果与前向散射能见度仪进行对比
%
% 殷振平
% 2024-01-21

clc;
close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\RawData\QingDao1';   % 数据主目录
tRange = [datenum(2023, 12, 25, 0, 0, 0), datenum(2023, 12, 31, 23, 59, 59)];   % 分析数据时间范围
distOffset = -48.75;   % 零点校正
iCh = 1;   % 读取通道下标（一共四个通道，其中有两个为近场和远场通道）
flagAve4Vis = false;   % 是否对每个扫描周期进行数据累加
iPrfInScan4Vis = 1;   % 如果不进行数据累加，则取每个扫描周期中的第iPrfInScan4Vis个廓线
flagOverlapCor = true;   % 是否进行重叠因子校正
olHeight = 1000;   % 盲区高度(米)，主要用于冼景洪算法
visRetMethod = 'xian';   % 能见度反演方法，建议使用冼景洪方法
overlapFile = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\wxzk\near-sea-shore-measurement\2023-12-29\overlap_20231229.mat';   % 重叠因子数据文件，重叠因子可以通过overlap_estimation_xxxx.m程序进行估算
visFile = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\wxzk\near-sea-shore-measurement\vis-fixed-platform.mat';   % 能见度数据文件

%% Read Data
fullData = struct();
fullData.mTime = [];
fullData.rawSignal = [];
fullData.nShots = [];
for thisDay = floor(tRange(1)):floor(tRange(2))
    dayPath = fullfile(dataFolder, datestr(thisDay, 'yyyy'), datestr(thisDay, 'mm'), datestr(thisDay, 'dd'));
    scanPaths = listdir(dayPath, '.*', 1);

    for iScan = 1:length(scanPaths)
        fprintf('Finished %6.2f%%: reading %s\n', (iScan - 1) / length(scanPaths) * 100, scanPaths{iScan});

        scanFiles = listfile(scanPaths{iScan}, '\w*', 1);
        if ~ flagAve4Vis
            scanFiles = scanFiles{iPrfInScan4Vis};
        end

        lData = readVIS(scanFiles, 'isDir', false);

        fullData.mTime = cat(2, fullData.mTime, mean(lData.startTime));
        fullData.nShots = cat(2, fullData.nShots, sum(lData.nShots));
        fullData.rawSignal = cat(1, fullData.rawSignal, sum(lData.rawSignal, 1));
    end
end

% Read Vis Data
visData = load(visFile);

%% Preprocess
range = ((1:lData.nBins(1)) - 0.5) * lData.hRes(1) + distOffset;

fullData.bg = nanmean(squeeze(fullData.rawSignal(:, iCh, (end - 20):end)), 2);
fullData.signal = squeeze(fullData.rawSignal(:, iCh, :)) - repmat(fullData.bg, 1, lData.nBins(1));
fullData.rcs = fullData.signal .* repmat(range, size(fullData.rawSignal, 1), 1).^2;
fullData.snr = fullData.signal ./ sqrt(squeeze(fullData.rawSignal(:, iCh, :)));
fullData.lowSNRMask = false(size(fullData.signal));
for iPrf = 1:size(fullData.signal, 1)
    snrPrf = fullData.snr(iPrf, :);
    snrPrf(range < olHeight) = NaN;
    snrLow = find(snrPrf < 1, 1);
    fullData.lowSNRMask(iPrf, snrLow:end) = true;
end

%% Overlap Correction
if flagOverlapCor
    ol = load(overlapFile);

    if exist('fullData', 'var')
        fullData.signal = fullData.signal ./ repmat(transpose(ol.ov), size(fullData.signal, 1), 1);
        fullData.rcs = fullData.rcs ./ repmat(transpose(ol.ov), size(fullData.rcs, 1), 1);
    end

    if exist('scanData', 'var')
        scanData.signal = scanData.signal ./ repmat(transpose(ol.ov), size(scanData.signal, 1), 1);
        scanData.rcs = scanData.rcs ./ repmat(transpose(ol.ov), size(scanData.rcs, 1), 1);
    end
end

%% extinction&visibility retrieval
fullData.ext = NaN(size(fullData.signal));
fullData.vis = NaN(size(fullData.signal));

for iPrf = 1:size(fullData.signal, 1)

    fullData.ext(iPrf, :) = extRet_Xian(range, fullData.signal(iPrf, :), fullData.bg(iPrf), 'minSNR', 0.5, 'rangeFullOverlap', 200);
end

fullData.ext(fullData.lowSNRMask) = NaN;
fullData.vis = ext2vis(fullData.ext);
fullData.vis(isnan(fullData.vis)) = 1e5;

%% Display

% full data rcs
figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
hold on;
rcsTmp = fullData.rcs;
rcsTmp(fullData.lowSNRMask) = NaN;
p1 = pcolor(fullData.mTime, range * 1e-3, transpose(rcsTmp));
p1.EdgeColor = 'none';
caxis([0, 0.1e11]);
colormap('jet');

xlabel('时间');
ylabel('距离 (千米)');
title(sprintf('激光测雾雷达时空信号图 (%s到%s)', datestr(fullData.mTime(1), 'yyyy-mm-dd'), datestr(fullData.mTime(end), 'mm-dd')));

xlim([fullData.mTime(1), fullData.mTime(end)]);
ylim([0, 4]);

set(gca, 'XMinorTick', 'on', 'xtick', floor(tRange(1)):ceil(tRange(2)), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 11);
ax = gca;
ax.XAxis.MinorTickValues = floor(tRange(1)):datenum(0, 1, 0, 6, 0, 0):ceil(tRange(2));
datetick(gca, 'x', 'mm-dd', 'keepticks', 'keeplimits');

cb = colorbar('Position', [0.89, 0.2, 0.03, 0.6], 'Units', 'Normalized');
titleHandle = get(cb, 'Title');
set(titleHandle, 'String', '[Mm-1sr-1]');
set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

% full data ext
figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
hold on;
extTmp = fullData.ext;
extTmp(fullData.lowSNRMask) = NaN;
p1 = pcolor(fullData.mTime, range * 1e-3, transpose(extTmp) * 1e3);
p1.EdgeColor = 'none';
caxis([0, 3]);
colormap('jet');

xlabel('时间');
ylabel('距离 (千米)');
title(sprintf('消光系数 (%s到%s)', datestr(fullData.mTime(1), 'yyyy-mm-dd'), datestr(fullData.mTime(end), 'mm-dd')));

xlim([fullData.mTime(1), fullData.mTime(end)]);
ylim([0, 4]);

set(gca, 'XMinorTick', 'on', 'xtick', floor(tRange(1)):ceil(tRange(2)), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 11);
ax = gca;
ax.XAxis.MinorTickValues = floor(tRange(1)):datenum(0, 1, 0, 6, 0, 0):ceil(tRange(2));
datetick(gca, 'x', 'mm-dd', 'keepticks', 'keeplimits');

cb = colorbar('Position', [0.89, 0.2, 0.03, 0.6], 'Units', 'Normalized');
titleHandle = get(cb, 'Title');
set(titleHandle, 'String', '[km-1]');
set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

% full data vis
figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
hold on;
visTmp = fullData.vis;
visTmp(fullData.lowSNRMask) = NaN;
p1 = pcolor(fullData.mTime, range * 1e-3, transpose(visTmp) * 1e-3);
p1.EdgeColor = 'none';
caxis([0, 30]);
load('vis_colormap.mat');
colormap(double(visColorbar) / 255);

xlabel('时间');
ylabel('距离 (千米)');
title(sprintf('能见度 (%s到%s)', datestr(fullData.mTime(1), 'yyyy-mm-dd'), datestr(fullData.mTime(end), 'mm-dd')));

xlim([fullData.mTime(1), fullData.mTime(end)]);
ylim([0, 4]);

set(gca, 'XMinorTick', 'on', 'xtick', floor(tRange(1)):ceil(tRange(2)), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 11);
ax = gca;
ax.XAxis.MinorTickValues = floor(tRange(1)):datenum(0, 1, 0, 6, 0, 0):ceil(tRange(2));
datetick(gca, 'x', 'mm-dd', 'keepticks', 'keeplimits');

cb = colorbar('Position', [0.89, 0.2, 0.03, 0.6], 'Units', 'Normalized');
titleHandle = get(cb, 'Title');
set(titleHandle, 'String', '[km]');
set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

% full data comparison
figure('Position', [0, 10, 550, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
s1 = plot(visData.mTime, visData.vis * 1e-3, 'color', 'k', 'Marker', '.', 'MarkerFaceColor', 'k', 'markeredgecolor', 'k', 'DisplayName', '能见度仪');
s2 = plot(fullData.mTime, fullData.vis(:, 80) * 1e-3, 'color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b', 'markeredgecolor', 'b', 'DisplayName', sprintf('激光测雾雷达(%3d米)', floor(range(80))));
hold off;

xlim([fullData.mTime(1), fullData.mTime(end)]);
ylim([0, 50]);

xlabel('时间');
ylabel('能见度 (千米)');
title(sprintf('能见度序列对比(%s到%s)', datestr(fullData.mTime(1), 'yyyy-mm-dd'), datestr(fullData.mTime(end), 'mm-dd')));

set(gca, 'XMinorTick', 'on', 'xtick', floor(tRange(1)):ceil(tRange(2)), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'fontsize', 11);
ax = gca;
ax.XAxis.MinorTickValues = floor(tRange(1)):datenum(0, 1, 0, 6, 0, 0):ceil(tRange(2));
datetick(gca, 'x', 'mm-dd', 'keepticks', 'keeplimits');
legend([s1, s2], 'location', 'NorthEast');