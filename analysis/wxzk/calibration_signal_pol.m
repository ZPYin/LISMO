% fog lidar with polarization

clc;
close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\12月26无人艇出海数据';   % input a scanning period. Path of scanning data.
tRange = [datenum(2023, 12, 26, 12, 0, 0), datenum(2023, 12, 26, 12, 20, 0)];   % temporal range for data input
distOffset = -52.5;   % distance offset. (m)
olHeight = 300;   % 
refDist = [3000, 3500];
calDist = [2000, 3000];
lr = 50;
maxRange = 10;   % km
overlapCor = false;
overlapFile = '';

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

lData = readVIS(dataFiles);

%% Preprocessing
range = ((1:lData.nBins(1)) - 0.5) * lData.hRes(1) + distOffset;
bg = nanmean(squeeze(lData.rawSignal(:, 3, (end - 100):(end - 20))), 2);
signal = squeeze(lData.rawSignal(:, 3, :)) - repmat(bg, 1, lData.nBins(1));
rcs = signal .* repmat(range, length(lData.hRes), 1).^2;
snr = signal ./ sqrt(squeeze(lData.rawSignal(:, 3, :)));

%% Interpolation of range corrected signal

% Display range corrected signal
figure('Position', [0, 30, 500, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = pcolor(lData.startTime, range * 1e-3, transpose(rcs));
p1.EdgeColor = 'none';
caxis([0, 2e9]);
colormap('jet');

xlabel('时间');
ylabel('距离 (千米)');
title(sprintf('激光测雾雷达时空信号图 (%s)', datestr(lData.startTime(1), 'yyyy-mm-dd')));

xlim(tRange);
ylim([0, 10]);

set(gca, 'XMinorTick', 'on', 'xtick', linspace(tRange(1), tRange(2), 6), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 10);
datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

colorbar();

%% Overlap Correction
meanRCS = nanmean(rcs, 1);
height = range .* cos(lData.zenithAng(1) / 180 * pi);
[mBsc, mExt] = MolModel(height, 1064, 'meteor', 'standard_atmosphere');   % Rayleigh Scattering
if overlapCor
end
mAttn = mBsc .* exp(-2 * nancumsum(mExt .* [range(1), diff(range)]));
isInRefDist = (range >= refDist(1)) & (range <= refDist(2));
ratioL2M = nansum(mAttn(isInRefDist)) / nansum(meanRCS(isInRefDist));

%% Retrieval

% calculate reference value
extSlopeMethod = movingslope(log(smooth(meanRCS, 8)), 8) / (-2) / (range(2) - range(1)) - transpose(mExt);
extRef = nanmean(extSlopeMethod(isInRefDist));
extRefStd = nanstd(extSlopeMethod(isInRefDist));

aBsc = fernald(range, nanmean(signal, 1), nansum(bg), lr, refDist, extRef / lr, mBsc, 4);

% display Rayleigh Fit
figure('Position', [0, 20, 500, 500], 'Units', 'Pixels', 'Color', 'w');

subfig = subfigPos([0.1, 0.1, 0.87, 0.8], 2, 1, 0, 0.1);

subplot('Position', subfig(1, :), 'Units', 'normalized');

hold on;

p1 = plot(range * 1e-3, mAttn, 'Color', [242, 89, 75] / 255, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', '分子信号');
p2 = plot(range * 1e-3, meanRCS * ratioL2M, 'Color', [65, 54, 89] / 255, 'LineWidth', 2, 'DisplayName', '雷达');

plot([refDist(1), refDist(1)] * 1e-3, [1e-20, 1e20], '--k');
plot([refDist(2), refDist(2)] * 1e-3, [1e-20, 1e20], '--k');

xlabel('水平距离 (千米)');
ylabel('衰减后向散射系数 (m-1sr-1)');

xlim([0, 5]);
ylim([1e-9, max(meanRCS * ratioL2M) * 4]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log');

legend([p1, p2], 'Location', 'SouthEast');

subplot('Position', subfig(2, :), 'Units', 'normalized');

hold on;

p1 = plot(range * 1e-3, mBsc * 1e6, 'Color', [242, 89, 75] / 255, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', '分子信号');
p2 = plot(range * 1e-3, aBsc * 1e6 * lr, 'Color', [65, 54, 89] / 255, 'LineWidth', 2, 'DisplayName', '雷达');

plot([refDist(1), refDist(1)] * 1e-3, [1e-20, 1e20], '--k');
plot([refDist(2), refDist(2)] * 1e-3, [1e-20, 1e20], '--k');

xlabel('水平距离 (千米)');
ylabel('消光系数 (m-1)');

xlim([0, 5]);
ylim([-0.5, 400]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

%% Lidar Calibration
aBscCor = aBsc;
aBscCor(range <= olHeight) = aBscCor(find(range >= olHeight, 1));
lc = meanRCS ./ ((aBscCor + mBsc) .* exp(-2 * nancumsum((aBscCor .* lr + mBsc) .* [range(1), diff(range)])));
isInCaliRange = (range >= calDist(1)) & (range <= calDist(2));
lcMean = nanmean(lc(isInCaliRange));
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

xlim([0, 5]);
ylim([0., lcMean * 3]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');