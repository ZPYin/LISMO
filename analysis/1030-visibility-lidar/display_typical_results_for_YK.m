% 制作用于意科公司能见度光量子雷达宣传使用的典型结果展示图
% 
% 作者：殷振平
% 日期：2026年05月27日
% 邮箱：zp.yin@whu.edu.cn

clc;
close all;
parentPath = fileparts(mfilename('fullpath'));

%% 参数定义
visSensorFile = 'vis-sensor-data.mat';   % 前向散射能见度仪数据文件
l0Folder = 'G:\backup\vis-lidar\20251221-wuhan\L0';   % 雷达原始数据文件目录
l1Folder = 'G:\backup\vis-lidar\20251221-wuhan\L1';   % 雷达产品文件目录
overlapFile = 'VIS_overlap_factor.txt';   % 重叠因子文件，重叠因子通过水平方法计算
tRange = [datenum(2025, 12, 22, 12, 0, 0), datenum(2025, 12, 24, 12, 0, 0)];
flagOLCor = true;   % 是否进行重叠因子修正
visDiffThresh = 2e4;   % [m]
dist_vis_lidar = 0.6;   % 前向散射能见度仪与激光雷达距离（千米）
flagDisplay = 'off';

%% Load Data

% vis sensor
load(visSensorFile);

% read overlap factor
olVal = [];
olHeight = [];
if exist(overlapFile, 'file') == 2
    fid = fopen(fullfile(parentPath, overlapFile), 'r');

    tmp = textscan(fid, '%f%f', 'HeaderLines', 1, 'Delimiter', ' ');
    olHeight = tmp{1};
    olVal = tmp{2};

    fclose(fid);
end

%% find dates
lidarSigMat = [];
mTimeMat = [];
visMat = [];

for thisDate = floor(tRange(1)):floor(tRange(2))

    % load level 0 data
    l0DataFile = fullfile(l0Folder, sprintf('%s_vis_lidar_l0.mat', datestr(thisDate, 'yyyy-mm-dd')));
    load(l0DataFile);

    % load level 1 data
    l1DataFile = fullfile(l1Folder, sprintf('%s_vis_lidar_l1_exp.mat', datestr(thisDate, 'yyyy-mm-dd')));
    load(l1DataFile);

    lidarSigMat = cat(2, lidarSigMat, lidarSig);
    mTimeMat = cat(2, mTimeMat, mTime);
    visMat = cat(2, visMat, visMat_Xian);

end

%% signal preprocessing
isInTRange = (mTimeMat >= tRange(1)) & (mTimeMat <= tRange(2));
lidarSigMat = lidarSigMat(:, isInTRange);
mTimeMat = mTimeMat(isInTRange);
thisBG = nanmean(lidarSigMat((end - 200):end, :), 1);
lidarSigCor = lidarSigMat - repmat(thisBG, size(thisBG, 1), 1);
if flagOLCor
    ol = interp1(olHeight, olVal, range);
    lidarSigCor = lidarSigCor ./ repmat(ol, 1, size(lidarSigCor, 2));
end
lidarRCS = lidarSigCor .* repmat(range, 1, size(lidarSigCor, 2)).^2;
visMat = visMat(:, isInTRange);
visMat(range < 400, :) = NaN;   % 质控，去除重叠因子影响区域内的结果

refIdx = find(range * 1e-3 >= dist_vis_lidar, 1);
visTimeseries = visMat(refIdx, :);
visDataInterp = interp1(visData.mTime, visData.vis, mTimeMat);

%% Display
figure('Position', [0, 30, 500, 450], 'Units', 'Pixels', 'Color', 'w', 'visible', flagDisplay);

% visibility comparison
subplot('Position', [0.1, 0.75, 0.77, 0.23], 'Units', 'Normalized');

hold on;
p2 = scatter(visData.mTime, visData.vis * 1e-3, 5, 'Marker', 's', 'MarkerEdgeColor', 'cyan', 'MarkerFaceColor', 'cyan', 'DisplayName', '前向散射能见度仪');
p1 = plot(mTimeMat, visTimeseries * 1e-3, '-r', 'DisplayName', '激光雷达');
hold off;

xlim(tRange);
ylim([0, 10]);

ylabel('能见度 (千米)');

text(0.03, 0.8, '(a) 能见度序列比对', 'Units', 'Normalized', 'FontSize', 13, 'Color', 'k', 'FontWeight', 'bold');

set(gca, 'YMinorTick', 'on', 'Box', 'on', 'XTick', linspace(min(mTimeMat), max(mTimeMat), 5), 'XTickLabel', '', 'layer', 'Top', 'LineWidth', 1, 'TickDir', 'out', 'FontSize', 11);

l = legend([p1, p2], 'Location', 'NorthEast');
l.FontSize = 10;

% range corrected signal
subplot('Position', [0.1, 0.44, 0.77, 0.27], 'Units', 'Normalized');

hold on;
p1 = pcolor(mTimeMat(1:10:end), range(1:4:end) * 1e-3, lidarRCS(1:4:end, 1:10:end) * 0.3e-10);
p1.EdgeColor = 'none';
hold off;

xlabel('');
ylabel('距离 (千米)');

xlim(tRange);
ylim([0, 10]);
caxis([0, 10]);
chilJetCmap = myColormap('chiljet');
colormap(gca, chilJetCmap);

text(0.03, 0.8, '(b) 能见度雷达信号', 'Units', 'Normalized', 'FontSize', 13, 'Color', 'k', 'FontWeight', 'bold');

set(gca, 'YMinorTick', 'on', 'Box', 'on', 'XTick', linspace(min(mTimeMat), max(mTimeMat), 5), 'XTickLabel', '', 'layer', 'Top', 'LineWidth', 1, 'TickDir', 'out', 'FontSize', 11);

cb = colorbar('Position', [0.92, 0.48, 0.03, 0.20], 'Units', 'Normalized');
titleHandle = get(cb, 'Title');
set(titleHandle, 'String', '[Mm^{-1}sr^{-1}]');
set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02, 'FontSize', 11);

% visibility colorplot
subplot('Position', [0.1, 0.12, 0.77, 0.28], 'Units', 'Normalized');

hold on;
p1 = pcolor(mTimeMat(1:10:end), range(1:4:end) * 1e-3, visMat(1:4:end, 1:10:end) * 1e-3);
p1.EdgeColor = 'none';

xlabel('时间 (日 时:分)');
ylabel('距离 (千米)');

xlim(tRange);
ylim([0, 10]);
caxis([0, 10]);
load('vis_colormap.mat');
colormap(gca, double(visColorbar) / 255);

text(0.03, 0.8, '(c) 能见度时空分布', 'Units', 'Normalized', 'FontSize', 13, 'Color', 'k', 'FontWeight', 'bold');
text(-0.05, -0.3, '2025年12月', 'Units', 'Normalized', 'FontSize', 13, 'Color', 'k');

set(gca, 'YMinorTick', 'on', 'Box', 'on', 'XTick', linspace(min(mTimeMat), max(mTimeMat), 5), 'layer', 'Top', 'LineWidth', 1, 'TickDir', 'out', 'FontSize', 11);

datetick(gca, 'x', 'dd HH:MM', 'keeplimits', 'keepticks');

cb = colorbar('Position', [0.92, 0.16, 0.03, 0.20], 'Units', 'Normalized');
titleHandle = get(cb, 'Title');
set(titleHandle, 'String', '[千米]');
set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02, 'FontSize', 11);

export_fig(gcf, fullfile(parentPath, 'whu_vis_lidar_demonstration_figure.png'), '-r300');
export_fig(gcf, fullfile(parentPath, 'whu_vis_lidar_demonstration_figure.pdf'), '-painters');
