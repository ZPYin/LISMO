%% 标定方法1：晴空标定

%% Parameter Definition
lidarPath = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\data\1030\2023-11-07';
tRnage = [datenum(2023, 11, 7, 13, 17, 0), datenum(2023, 11, 7, 13, 18, 0)];
hOffset = -285;
olHeight = 1000;   % 
deadTime = [20, 20];   % 死时间，单位ns
refDist = [12000, 16000];
calDist = [2000, 4000];
lr = 50;
maxRange = 10;   % km
overlapCor = false;
overlapFile = '';

%% Read Data
lData = readALADats(lidarPath, 'tRange', tRange, 'nmaxbin', 2100);

%% Preprocessing

% 死时间校正
rawSig = lData.rawSignal;
pc2cr = 1 / (1500 / lData.hRes(1) * lData.nShots(1));
rawSigCR = rawSig * pc2cr;
rawSigCRCor = rawSigCR ./ (1 - repmat(reshape(transpose(deadTime), [], 1, 1), 1, size(rawSigCR, 2), size(rawSigCR, 3)) * 1e-9 .* rawSigCR);
rawSigPCCor = rawSigCRCor / pc2cr;
range = (1:size(rawSigPCCor, 2)) * lData.hRes(1) + hOffset;

% 背景扣除
bg = nanmean(rawSigCRCor(:, (end - 10):end, :), 2);
sigNoBg = rawSigPCCor - repmat(bg, 1, size(rawSigCRCor, 2), 1);
rcs = squeeze(sigNoBg(2, :, :)) .* repmat(range', 1, length(lData.mTime)).^2;

%% Overlap Correction
meanRCS = nanmean(rcs, 2);
height = range .* cos(88 / 180 * pi);
[mBsc, mExt] = MolModel(height, 1030, 'meteor', 'standard_atmosphere');   % Rayleigh Scattering
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

aBsc = fernald(range, nanmean(squeeze(sigNoBg(2, :, :)), 2), nansum(bg(2, 1, :)), lr, refDist, extRef / lr, mBsc, 4);

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

xlim([0, 20]);
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

xlim([0, 20]);
ylim([-0.5, 150]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

%% Lidar Calibration
aBscCor = aBsc;
aBscCor(range <= olHeight) = aBscCor(find(range >= olHeight, 1));
lc = meanRCS' ./ ((aBscCor + mBsc) .* exp(-2 * nancumsum((aBscCor .* lr + mBsc) .* [range(1), diff(range)])));
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

xlim([0, 20]);
ylim([0., lcMean * 3]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

%% 标定方法2：前散标定