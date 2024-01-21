% fog lidar with polarization

clc;
% close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\海试数据、效果图\10月21长航实验原始数据';
tRange = [datenum(2023, 10, 21, 6, 30, 0), datenum(2023, 10, 21, 21, 14, 0)];
tCaliRange = [datenum(2023, 10, 21, 18, 30, 0), datenum(2023, 10, 21, 18, 32, 0)];
distOffset = 55.75;
olHeight = 300;
refDist = [1500, 3000];
calDist = [1500, 2000];
lr = 30;
maxRange = 10;   % km
fitRange = [0, 4];   % signal glue range. (MHz)
visRetMethod = 'xian';
flagReadData = false;
overlapCor = false;
overlapFile = 'overlap_20231021.mat';
visFile = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\wxzk\shipborne-measurements\visibility.mat';

%% Read Data
if flagReadData
    lData = readVIS(dataFolder, 'tRange', tRange, 'isDir', true, 'debug', true);
end

% read visibility
visData = load(visFile);
load('vis_colormap.mat');

%% Preprocessing
range = ((1:lData.nBins(1)) - 0.5) * lData.hRes(1) - distOffset;
sigFR = squeeze(lData.rawSignal(:, 3, :));
sigNR = squeeze(lData.rawSignal(:, 4, :));
sigFRNoBG = sigFR - nanmean(sigFR(:, (end - 15):end), 2);
sigNRNoBG = sigNR - nanmean(sigNR(:, (end - 15):end), 2);
% isSaturated = (sigFR / (lData.hRes(1) / 15 * lData.nShots(1) * 100 * 1e-3) > fitRange(2));
isSaturated = false(size(sigFR));
signal = sigFR;
signal(isSaturated) = ((sigNR(isSaturated) / (lData.hRes(1) / 15 * lData.nShots(1) * 100 * 1e-3)) * 52.9573 - 0.0060) * (lData.hRes(1) / 15 *lData.nShots(1) * 100 * 1e-3);
bg = nanmean(signal(:, (end - 10):end), 2);
signal = signal - repmat(bg, 1, lData.nBins(1));
rcs = signal .* repmat(range, length(lData.hRes), 1).^2;
snr = signal ./ sqrt(signal + repmat(bg, 1, lData.nBins(1)));
lowSNRMask = false(size(signal));
for iPrf = 1:size(rcs, 1)
    snrPrf = snr(iPrf, :);
    snrPrf(range < olHeight) = NaN;
    snrLow = find(snrPrf < 1, 1);
    lowSNRMask(iPrf, snrLow:end) = true;
end

%% Overlap Correction
if overlapCor
    ol = load(overlapFile);
    rcs = rcs ./ repmat(transpose(ol.ov), size(rcs, 1), 1);
end

lcMean = 0.3e15;

%% visibility
ext = NaN(size(signal));
for iPrf = 1:length(lData.hRes)

    if lData.startTime(iPrf) <= datenum(2023, 10, 21, 20, 30, 0)
        ext(iPrf, :) = extRet_Xian(range, signal(iPrf, :), bg(iPrf), 'minSNR', 1, 'rangeFullOverlap', 800);
    else
        [~, ext(iPrf, :)] = extRet_Holger(range, signal(iPrf, :), ...
            'calibration_constant', lcMean, ...
            'fullOverlapR', 300, ...
            'elevation_angle', lData.zenithAng(iPrf));
    end

end

ext(lowSNRMask) = NaN;
vis = ext2vis(ext);
vis(vis > 1e5) = 1e5;

%% Display Vis

% far-range
figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
hold on;
rcsTmp = rcs;
rcsTmp(lowSNRMask) = NaN;
p1 = pcolor(lData.startTime, range * 1e-3, transpose(vis) * 1e-3);
p1.EdgeColor = 'none';
caxis([0, 30]);
colormap(double(visColorbar) / 255);

xlabel('时间');
ylabel('距离 (千米)');
title(sprintf('能见度时空分布图 (%s)', datestr(lData.startTime(1), 'yyyy-mm-dd')));

xlim(tRange);
ylim([0, 10]);

set(gca, 'XMinorTick', 'on', 'xtick', linspace(tRange(1), tRange(2), 6), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 11);
datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

cb = colorbar('Position', [0.89, 0.2, 0.03, 0.6], 'Units', 'Normalized');
titleHandle = get(cb, 'Title');
set(titleHandle, 'String', '[千米]');
set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

%% Display visiblity
figure('Position', [0, 10, 550, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
s1 = plot(visData.mTime, visData.vis * 1e-3, 'color', 'k', 'Marker', '.', 'MarkerFaceColor', 'k', 'markeredgecolor', 'k', 'DisplayName', '能见度仪');
s2 = plot(lData.startTime, vis(:, 70) * 1e-3, 'color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b', 'markeredgecolor', 'b', 'DisplayName', sprintf('激光测雾雷达(%3d米)', floor(range(70))));
hold off;

xlim(tRange);
ylim([0, 50]);

xlabel('时间');
ylabel('能见度 (千米)');
title(sprintf('能见度序列对比(%s)', datestr(lData.startTime(1), 'yyyy-mm-dd')));

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'xtick', linspace(tRange(1), tRange(2), 6), 'Box', 'on', 'FontSize', 11);
datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');
legend([s1, s2], 'location', 'NorthEast');