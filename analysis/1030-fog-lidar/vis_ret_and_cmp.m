clc;
close all;
global LISMO_VARS;

%% 参数定义
lidarPath = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\data\1030\2023-11-07';
tRange = [datenum(2023, 11, 7, 11, 40, 0), datenum(2023, 11, 7, 12, 50, 0)];
deadTime = [20, 20];   % 死时间，单位ns
hOffset = -225;
zenithAngle = 0;
dbFile = '';
cloudMinExt = 200e-6;   % m-1
overlapFile = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\data\1030-fog-lidar-overlap-2023-11-07.mat';
idxVisCmp = 387;   % 用于对比的能见度距离门下标
ratioNear = 76;

%% 数据读取
lData = readALADats(lidarPath, 'tRange', tRange, 'nmaxbin', 2000);
ol = load(overlapFile);

%% 数据预处理

% 死时间校正
rawSig = lData.rawSignal;
pc2cr = 1 / (1500 / lData.hRes(1) * lData.nShots(1));
rawSigCR = rawSig * pc2cr;
rawSigCRCor = rawSigCR ./ (1 - repmat(reshape(transpose(deadTime), [], 1, 1), 1, size(rawSigCR, 2), size(rawSigCR, 3)) * 1e-9 .* rawSigCR);
rawSigPCCor = rawSigCRCor / pc2cr;
range = (1:size(rawSigPCCor, 2)) * lData.hRes(1) + hOffset;

% 背景扣除
bg = nanmean(rawSigCRCor(:, (end - 50):end, :), 2);
sigNoBg = rawSigPCCor - repmat(bg, 1, size(rawSigCRCor, 2), 1);

% 重叠因子修正
sigNearNoBg = squeeze(sigNoBg(1, :, :));
sigFarNoBg = squeeze(sigNoBg(2, :, :));
overlapNear = interp1(ol.range, ol.overlapNear, range);
overlapFar = interp1(ol.range, ol.overlapFar, range);
sigFarNoBgCor = overlapCor(range, sigFarNoBg, overlapFar, 'glueRange', [1000, 1400]);
sigNearNoBgCor = sigNearNoBg ./ repmat(transpose(overlapNear), 1, size(sigNearNoBg, 2)) * ratioNear;

% 信号拼接
sigMerge = signalMerge(sigFarNoBgCor, sigNearNoBgCor, range, [0, 100], 1, 0);

% 信噪比
snrNear = sigNearNoBg ./ sqrt(squeeze(rawSigPCCor(1, :, :)));
snrFar = sigFarNoBg ./ sqrt(squeeze(rawSigPCCor(2, :, :)));
isNoisy = (snrNear < 1) & (snrFar < 1);

%% 雷达标定
lc = getLC(lData.mTime, dbFile);
attnBsc = (sigMerge .* repmat(reshape(range, size(sigMerge, 1), 1), 1, size(sigNearNoBg, 2)).^2) / lc;
attnBsc(range <= 20, :) = NaN;

%% 消光系数反演（前向反演）
fExt = zeros(size(attnBsc));
feat = zeros(size(attnBsc));   % 0: no-cloud; 1: cloud; 2: unknown
lr = ones(size(attnBsc)) * 50;

[mBsc, mExt] = MolModel(range .* cos(zenithAngle / 180 * pi), 1030, 'meteor', 'standard_atmosphere');

isConverge = false;
counter = 0;
while ((~ isConverge) && (counter < 10))
    % 前向反演
    [fBsc, fExt] = quasiRetrieval(range, attnBsc, repmat(mExt', 1, size(lr, 2)), repmat(mBsc', 1, size(lr, 2)), lr, 'nIters', 10, 'flagAutoConverge', true);

    % 特征识别
    featBefore = feat;
    isCloud = (fExt >= cloudMinExt);
    isAerosol = (fExt < cloudMinExt);
    feat(isCloud) = 1;
    feat(isAerosol) = 0;
    feat(isNoisy) = 2;
    lr(isCloud) = 20;
    lr(isAerosol) = 45;
    lr(isNoisy) = 45;

    isSameFeatFrac = sum(sum(feat == featBefore, 2)) / numel(feat);
    isConverge = (isSameFeatFrac > 0.99);
    counter = counter + 1;
end

%% 能见度反演
vis = ext2vis((fExt + repmat(mExt', 1, size(lr, 2))) * (1030 / 532) .^ 1.7);

% 数据可视化
figure('Position', [0, 20, 500, 600], 'Units', 'Pixels', 'Color', 'w');

subfig = subfigPos([0.1, 0.1, 0.83, 0.86], 3, 1, 0, 0.05);

% 距离修正信号
subplot('Position', subfig(1, :), 'Units', 'normalized');
hold on;
p1 = pcolor(lData.mTime, range * 1e-3, attnBsc);
p1.EdgeColor = 'none';

xlabel('');
ylabel('距离 (千米)');
title('距离修正信号');

xlim(tRange);
ylim([0, 10]);
caxis([0, 2e-6]);
colormap('jet');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'layer', 'top', 'Tickdir', 'out', 'XTickLabel', '', 'XTick', linspace(tRange(1), tRange(2), 5));

colorbar();

% 消光系数
subplot('Position', subfig(2, :), 'Units', 'normalized');
hold on;
p1 = pcolor(lData.mTime, range * 1e-3, fExt);
p1.EdgeColor = 'none';

xlabel('');
ylabel('距离 (千米)');
title('消光系数 (m-1)');

xlim(tRange);
ylim([0, 10]);
caxis([0, 1e-4]);
colormap('jet');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'layer', 'top', 'Tickdir', 'out', 'XTickLabel', '', 'XTick', linspace(tRange(1), tRange(2), 5));

colorbar();

% 能见度
subplot('Position', subfig(3, :), 'Units', 'normalized');
hold on;
p1 = pcolor(lData.mTime, range * 1e-3, vis);
p1.EdgeColor = 'none';

xlabel('时间');
ylabel('距离 (千米)');
title('能见度 (米)');

xlim(tRange);
ylim([0, 10]);
caxis([0, 5e4]);
colormap(gca, 'hot');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'layer', 'top', 'Tickdir', 'out', 'XTick', linspace(tRange(1), tRange(2), 5));
datetick(gca, 'x', 'HH:MM', 'Keeplimits', 'keepticks');

colorbar();

%% 能见度对比

% 读取前向散射能见度仪数据
visDataFile = fullfile(LISMO_VARS.projectDir, 'data', 'vis35.mat');
a = load(visDataFile, 'vis35');
vis35 = a.vis35;
isInTRange = (vis35.mTime >= (tRange(1) + datenum(0, 1, 0, 8, 0, 0))) & (vis35.mTime <= (tRange(2) + datenum(0, 1, 0, 8, 0, 0)));
vis35Time = vis35.mTime(isInTRange);
vis1min = vis35.vis1min(isInTRange);
vis1minInterp = interp1(vis35Time, vis1min, lData.mTime + datenum(0, 1, 0, 8, 0, 0));
vis1minRelBias = abs(vis1minInterp - vis(idxVisCmp, :)) ./ vis1minInterp;
meanVisBias = nanmean(vis1minRelBias);

figure('Position', [0, 20, 600, 400], 'Units', 'Pixels', 'Color', 'w');

subfig = subfigPos([0.1, 0.1, 0.83, 0.8], 2, 1, 0, 0.1);

% 能见度时间序列
subplot('Position', subfig(1, :), 'Units', 'normalized');

hold on;

p1 = plot(vis35Time, vis1min, 'Marker', '.', 'LineStyle', 'none', 'Color', 'k', 'DisplayName', '前向散射能见度仪');
p2 = plot(lData.mTime + datenum(0, 1, 0, 8, 0, 0), vis(idxVisCmp, :), 'Marker', 's', 'LineStyle', 'none', 'Color', 'b', 'DisplayName', '激光测雾雷达');

xlabel('');
ylabel('能见度 (米)');
title('能见度序列对比');

xlim(tRange + datenum(0, 1, 0, 8, 0, 0));
ylim([0, 50000]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');
datetick(gca, 'x', 'HH:MM', 'Keepticks', 'Keeplimits');
legend([p1, p2], 'Location', 'SouthEast');


% 能见度偏差
subplot('Position', subfig(2, :), 'Units', 'normalized');

hold on;
p1 = plot(lData.mTime + datenum(0, 1, 0, 8, 0, 0), vis1minRelBias * 100, 'Color', 'b', 'DisplayName', '激光测雾雷达');
p2 = plot(lData.mTime + datenum(0, 1, 0, 8, 0, 0), ones(size(vis1minRelBias)) * 30, '--k');
p3 = plot(lData.mTime + datenum(0, 1, 0, 8, 0, 0), ones(size(vis1minRelBias)) * meanVisBias * 100, '--b', 'LineWidth', 2);
hold off;

text(0.6, 0.8, sprintf('平均相对偏差为%5.2f%%', meanVisBias * 100), 'Units', 'Normalized');

xlabel('');
ylabel('相对偏差 (%)');

xlim(tRange + datenum(0, 1, 0, 8, 0, 0));
ylim([0, 50]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');
datetick(gca, 'x', 'HH:MM', 'Keepticks', 'Keeplimits');