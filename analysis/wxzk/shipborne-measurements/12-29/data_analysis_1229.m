% fog lidar with polarization

clc;
close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\海试数据、效果图\无人艇12月18、29数据\12.29';   % input a scanning period. Path of scanning data.
tRange = [datenum(2023, 12, 29, 8, 46, 0), datenum(2023, 12, 29, 16, 43, 0)];   % temporal range for data input
tCaliRange = [datenum(2023, 12, 29, 11, 0, 0), datenum(2023, 12, 29, 11, 10, 0)];
distOffset = 55.75;
olHeight = 300;   % 
refDist = [800, 1000];
calDist = [300, 400];
lr = 30;
maxRange = 10;   % km
visRetMethod = 'xian';
flagReadData = false;
overlapCor = true;
overlapFile = 'overlap_20231228.mat';
visFile = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\wxzk\shipborne-measurements\visibility.mat';

%% Read Data
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

if flagReadData
    lData = readVIS(dataFiles, 'debug', true);
end

% read visibility
visData = load(visFile);

%% Preprocessing

%% Preprocessing
startTime = lData.startTime(1):datenum(0, 1, 0, 0, 1, 0):lData.startTime(end);
range = ((1:lData.nBins(1)) - 0.5) * lData.hRes(1) - distOffset;
bg = NaN(1, length(startTime));
signal = NaN(length(startTime), length(range));

for iPrf = 1:length(lData.startTime)
    [~, tIdx] = min(abs(startTime - lData.startTime(iPrf)));
    bg(tIdx) = nanmean(squeeze(lData.rawSignal(iPrf, 3, (end - 100):(end - 20))));
    signal(tIdx, :) = transpose(squeeze(lData.rawSignal(iPrf, 3, :))) - repmat(bg(tIdx), 1, lData.nBins(1));
end
rcs = signal .* repmat(range, length(startTime), 1).^2;
snr = signal ./ (signal + repmat(bg', 1, length(range)));
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

%% Interpolation of range corrected signal

% Display range corrected signal
figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
hold on;
rcsTmp = rcs;
rcsTmp(lowSNRMask) = NaN;
p1 = pcolor(startTime, range * 1e-3, transpose(rcsTmp));
p1.EdgeColor = 'none';
caxis([0, 3e9]);
colormap('jet');

xlabel('时间');
ylabel('距离 (千米)');
title(sprintf('激光测雾雷达时空信号图 (%s)', datestr(startTime(1), 'yyyy-mm-dd')));

xlim(tRange);
ylim([0, 10]);

set(gca, 'XMinorTick', 'on', 'xtick', linspace(tRange(1), tRange(2), 6), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 11);
datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

cb = colorbar('Position', [0.89, 0.2, 0.03, 0.6], 'Units', 'Normalized');
titleHandle = get(cb, 'Title');
set(titleHandle, 'String', '[a.u.]');
set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

%% lidar calibration
isInCaliRange = (startTime >= tCaliRange(1)) & (startTime <= tCaliRange(2));
meanRCS = nanmean(rcs(isInCaliRange, :), 1);
height = range .* cos(lData.zenithAng(1) / 180 * pi);
[mBsc, mExt] = MolModel(height, 1064, 'meteor', 'standard_atmosphere');   % Rayleigh Scattering
mAttn = mBsc .* exp(-2 * nancumsum(mExt .* [range(1), diff(range)]));
isInRefDist = (range >= refDist(1)) & (range <= refDist(2));
ratioL2M = nansum(mAttn(isInRefDist)) / nansum(meanRCS(isInRefDist));

%% Retrieval

% calculate reference value
extSlopeMethod = movingslope(log(smooth(meanRCS, 8)), 8) / (-2) / (range(2) - range(1)) - transpose(mExt);
extRef = nanmean(extSlopeMethod(isInRefDist));
extRefStd = nanstd(extSlopeMethod(isInRefDist));

aBsc = fernald(range, nanmean(signal(isInCaliRange, :), 1), nansum(bg(isInCaliRange)), lr, refDist, extRef / lr, mBsc, 4);

% display Rayleigh Fit
figure('Position', [0, 20, 500, 500], 'Units', 'Pixels', 'Color', 'w');

subfig = subfigPos([0.12, 0.1, 0.85, 0.8], 2, 1, 0, 0.1);

subplot('Position', subfig(1, :), 'Units', 'normalized');

hold on;

p1 = plot(range * 1e-3, mAttn * 1e3, 'Color', [242, 89, 75] / 255, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', '分子信号');
p2 = plot(range * 1e-3, meanRCS * ratioL2M * 1e3, 'Color', [65, 54, 89] / 255, 'LineWidth', 2, 'DisplayName', '雷达');

plot([refDist(1), refDist(1)] * 1e-3, [1e-20, 1e20], '--k');
plot([refDist(2), refDist(2)] * 1e-3, [1e-20, 1e20], '--k');

xlabel('水平距离 (千米)');
ylabel('衰减后向散射系数 (km-1sr-1)');

xlim([0, 5]);
ylim([1e-6, 1e-3]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log');

legend([p1, p2], 'Location', 'SouthEast');

subplot('Position', subfig(2, :), 'Units', 'normalized');

hold on;

p1 = plot(range * 1e-3, mBsc * 1e3, 'Color', [242, 89, 75] / 255, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', '分子信号');
p2 = plot(range * 1e-3, aBsc * 1e3 * lr, 'Color', [65, 54, 89] / 255, 'LineWidth', 2, 'DisplayName', '雷达');

plot([refDist(1), refDist(1)] * 1e-3, [1e-20, 1e20], '--k');
plot([refDist(2), refDist(2)] * 1e-3, [1e-20, 1e20], '--k');

xlabel('水平距离 (千米)');
ylabel('消光系数 (km-1)');

xlim([0, 5]);
ylim([-0.001, 1]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

%% Lidar Calibration
aBscCor = aBsc;
aBscCor(range <= olHeight) = aBscCor(find(range >= olHeight, 1));
lc = meanRCS ./ ((aBscCor + mBsc) .* exp(-2 * nancumsum((aBscCor .* lr + mBsc) .* [range(1), diff(range)])));
isInCaliRange = (range >= calDist(1)) & (range <= calDist(2));
lcMean = 1e14;
lcStd = nanstd(lc(isInCaliRange));
meanL2M = nanmean(meanRCS(isInCaliRange)) / nanmean(mAttn(isInCaliRange));

figure('Position', [0, 30, 500, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
plot(range * 1e-3, lc, '-k', 'LineWidth', 2);
plot([calDist(1), calDist(1)] * 1e-3, [1e-20, 1e20], '--k');
plot([calDist(2), calDist(2)] * 1e-3, [1e-20, 1e20], '--k');

text(0.3, 0.8, sprintf('雷达标定常数平均: %5.3e\\pm%5.3e\n参考分子标定常数:%5.3e\n', lcMean, lcStd, meanL2M), 'Units', 'normalized');

xlabel('水平距离 (千米)');
ylabel('雷达标定常数 (每单廓线)');

xlim([0, 3]);
ylim([0., lcMean * 3]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

% Display attenuated backscatter (far-range)
figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
hold on;
rcsTmp = rcs;
rcsTmp(lowSNRMask) = NaN;
p1 = pcolor(startTime, range * 1e-3, transpose(rcsTmp) / lcMean * 1e6);
p1.EdgeColor = 'none';
caxis([0, 20]);
colormap('jet');

xlabel('时间');
ylabel('距离 (千米)');
title(sprintf('激光测雾雷达时空信号图 (%s)', datestr(startTime(1), 'yyyy-mm-dd')));

xlim(tRange);
ylim([0, 10]);

set(gca, 'XMinorTick', 'on', 'xtick', linspace(tRange(1), tRange(2), 6), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 11);
datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

cb = colorbar('Position', [0.89, 0.2, 0.03, 0.6], 'Units', 'Normalized');
titleHandle = get(cb, 'Title');
set(titleHandle, 'String', '[Mm-1sr-1]');
set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

% Display attenuated backscatter (near-range)
figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
hold on;
rcsTmp = rcs;
rcsTmp(lowSNRMask) = NaN;
p1 = pcolor(startTime, range * 1e-3, transpose(rcsTmp) / lcMean * 1e6);
p1.EdgeColor = 'none';
caxis([0, 20]);
colormap('jet');

xlabel('时间');
ylabel('距离 (千米)');
title(sprintf('激光测雾雷达时空信号图 (%s)', datestr(startTime(1), 'yyyy-mm-dd')));

xlim(tRange);
ylim([0, 4]);

set(gca, 'XMinorTick', 'on', 'xtick', linspace(tRange(1), tRange(2), 6), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 11);
datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

cb = colorbar('Position', [0.89, 0.2, 0.03, 0.6], 'Units', 'Normalized');
titleHandle = get(cb, 'Title');
set(titleHandle, 'String', '[Mm-1sr-1]');
set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

%% visibility
ext = NaN(size(signal));
for iPrf = 1:length(lData.hRes)
    if isnan(bg(iPrf))
        continue;
    end

    if strcmpi(visRetMethod, 'xian')
        ext(iPrf, :) = extRet_Xian(range, signal(iPrf, :), bg(iPrf), 'minSNR', 0.1, 'rangeFullOverlap', 150);
    elseif strcmpi(visRetMethod, 'quasi')
        [~, ext(iPrf, :)] = extRet_Holger(range, signal(iPrf, :), ...
            'calibration_constant', lcMean, ...
            'fullOverlapR', 100, ...
            'elevation_angle', lData.zenithAng(iPrf));
    else
        warning('Undefined extinction retrieval method.');
    end

end

ext(lowSNRMask) = NaN;
vis = ext2vis(ext);
vis(isnan(vis)) = NaN;

%% Display Vis

% far-range
figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
hold on;
rcsTmp = rcs;
rcsTmp(lowSNRMask) = NaN;
p1 = pcolor(startTime, range * 1e-3, transpose(vis) * 1e-3);
p1.EdgeColor = 'none';
caxis([0, 30]);
colormap(double(visColorbar) / 255);

xlabel('时间');
ylabel('距离 (千米)');
title(sprintf('能见度时空分布图 (%s)', datestr(startTime(1), 'yyyy-mm-dd')));

xlim(tRange);
ylim([0, 10]);

set(gca, 'XMinorTick', 'on', 'xtick', linspace(tRange(1), tRange(2), 6), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 11);
datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

cb = colorbar('Position', [0.89, 0.2, 0.03, 0.6], 'Units', 'Normalized');
titleHandle = get(cb, 'Title');
set(titleHandle, 'String', '[千米]');
set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

% near-range
figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
hold on;
rcsTmp = rcs;
rcsTmp(lowSNRMask) = NaN;
p1 = pcolor(startTime, range * 1e-3, transpose(vis) * 1e-3);
p1.EdgeColor = 'none';
caxis([0, 30]);
colormap(double(visColorbar) / 255);

xlabel('时间');
ylabel('距离 (千米)');
title(sprintf('能见度时空分布图 (%s)', datestr(startTime(1), 'yyyy-mm-dd')));

xlim(tRange);
ylim([0, 4]);

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
s2 = plot(startTime, vis(:, 30) * 1e-3, 'color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b', 'markeredgecolor', 'b', 'DisplayName', sprintf('激光测雾雷达(%3d米)', floor(range(30))));
hold off;

xlim(tRange);
ylim([0, 50]);

xlabel('时间');
ylabel('能见度 (千米)');
title(sprintf('能见度序列对比(%s)', datestr(startTime(1), 'yyyy-mm-dd')));

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);
datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');
legend([s1, s2], 'location', 'NorthEast');