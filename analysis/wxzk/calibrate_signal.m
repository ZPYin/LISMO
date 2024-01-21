clc;
close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\RawData\QingDao\2023\12\27\H_SCAN15_135_2_20231227085349';   % input a scanning period. Path of scanning data.
tRnage = [datenum(2023, 12, 25, 0, 0, 0), datenum(2023, 12, 25, 0, 30, 0)];   % temporal range for data input
distOffset = -52.5;   % distance offset. (m)
olHeight = 300;   % 
refDist = [2200, 3000];
calDist = [1000, 2000];
lr = 50;
maxRange = 5;   % km
overlapCor = false;
overlapFile = '';

%% Read Data
dataFilesSearched = listfile(dataFolder, '\w*.VIS', 1);
dataFiles = cell(0);
for iFile = 1:length(dataFilesSearched)
    fileBasename = basename(dataFilesSearched{iFile});
    fileTime = datenum(fileBasename(3:(end - 4)), 'yyyymmddHHMMSS');

    isInTRange = (fileTime >= tRnage(1)) & (fileTime < tRnage(2));
    if isInTRange
        dataFiles = cat(2, dataFiles, dataFilesSearched{iFile});
    end
end

lData = readVIS(dataFilesSearched);

%% Preprocessing
range = ((1:lData.nBins(1)) - 0.5) * lData.hRes(1) + distOffset;
bg = nanmean(squeeze(lData.rawSignal(:, 1, (end - 100):(end - 20))), 2);
signal = squeeze(lData.rawSignal(:, 1, :)) - repmat(bg, 1, lData.nBins(1));
rcs = signal .* repmat(range, length(lData.hRes), 1).^2;
snr = signal ./ sqrt(squeeze(lData.rawSignal(:, 1, :)));

%% Interpolation of range corrected signal
x = transpose(range) * (cos((90 - lData.zenithAng) / 180 * pi) .* sin(lData.azimuthAng / 180 * pi));
y = transpose(range) * (cos((90 - lData.zenithAng) / 180 * pi) .* cos(lData.azimuthAng / 180 * pi));

% Display range corrected signal
figure('Position', [0, 30, 700, 300], 'Units', 'Pixels', 'Color', 'w');

subfig = subfigPos([0.1, 0.15, 0.8, 0.75], 1, 2, 0.1, 0);

subplot('Position', subfig(1, :), 'Units', 'normalized');

hold on;
p1 = pcolor(x * 1e-3, y * 1e-3, transpose(rcs));
p1.EdgeColor = 'none';
caxis([0, 6e9]);
colormap('jet');

xlabel('x轴距离 (千米)');
ylabel('y轴距离 (千米)');
title(sprintf('激光测雾雷达时空信号图 (%s)', datestr(lData.startTime(1), 'yyyy-mm-dd HH:MM')));

xlim([-maxRange, maxRange]);
ylim([-maxRange, maxRange]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'plotboxaspectratio', [1, 1, 1], 'fontsize', 10);

colorbar('Position', [0.45, 0.3, 0.02, 0.4]);

subplot('Position', subfig(2, :), 'Units', 'normalized');

hold on;
p1 = pcolor(x * 1e-3, y * 1e-3, transpose(snr));
p1.EdgeColor = 'none';
caxis([0, 10]);
colormap('jet');

xlabel('x轴距离 (千米)');
ylabel('y轴距离 (千米)');
title('信噪比');

xlim([-maxRange, maxRange]);
ylim([-maxRange, maxRange]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'plotboxaspectratio', [1, 1, 1], 'fontsize', 10);

colorbar('Position', [0.9, 0.3, 0.02, 0.4]);

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

xlim([0, maxRange]);
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

xlim([0, maxRange]);
ylim([-0.5, 500]);

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

xlim([0, maxRange]);
ylim([0., lcMean * 3]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');