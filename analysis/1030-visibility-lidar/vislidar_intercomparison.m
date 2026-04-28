% 用于两台能见度激光雷达数据对比
%
% 作者：殷振平
% 邮箱：zp.yin@whu.edu.cn
% 日期：2026-04-28

clc;
close all;

%% Parameter Definition
dataPath1 = 'G:\backup\vis-lidar\20260428-tianjin-evaluation\一致性数据\VIS1';   % 数据目录1
dataPath2 = 'G:\backup\vis-lidar\20260428-tianjin-evaluation\一致性数据\VIS2';   % 数据目录2
savePath = 'C:\Users\zhenp\OneDrive\Desktop';   % 图片结果输出目录
dataStr = '20260425';
flagOLCor = false;   % 是否进行重叠因子修正
hRange = [0, 10];
olFile = 'overlap_20250101.mat';   % 重叠因子文件
hFullOL = 500;   % 完全进入视场高度（m）
distOffset = -55;   % 预触发点个数（可通过信号第一个峰值进行判断）
visible = 'on';   % 是否进行结果可视化

%% Read Data
mDate = datenum(dataStr, 'yyyymmdd');
lData1 = readALADats(fullfile(dataPath1, datestr(mDate, 'yyyymmdd')), 'datatype', 2);
lData2 = readALADats(fullfile(dataPath2, datestr(mDate, 'yyyymmdd')), 'datatype', 2);

%% Preprocessing
range1 = transpose(((1:size(lData1.rawSignal, 2)) - 0.5 + distOffset) * lData1.hRes(1));
bg1 = mean(lData1.rawSignal(:, (end - 30):(end - 5), :), 2);
noise1 = std(lData1.rawSignal(:, (end - 30):(end - 5), :), 0, 2);
signal1 = squeeze(lData1.rawSignal - repmat(bg1, size(lData1.rawSignal, 1), size(lData1.rawSignal, 2), 1));
mTime1 = lData1.mTime;
if flagOLCor
    ol = load(olFile);
    signal = signal ./ repmat(transpose(ol.ov), 1, size(lData1.rawSignal, 2));
end
rcs1 = signal1 .* repmat(range1, 1, size(lData1.rawSignal, 3)).^2;   % 距离修正信号
snr1 = signal1 ./ repmat(transpose(squeeze(noise1)), size(lData1.rawSignal, 2), 1);   % 信噪比

range2 = transpose(((1:size(lData2.rawSignal, 2)) - 0.5 + distOffset) * lData2.hRes(1));
bg2 = mean(lData2.rawSignal(:, (end - 30):(end - 5), :), 2);
noise2 = std(lData2.rawSignal(:, (end - 30):(end - 5), :), 0, 2);
signal2 = squeeze(lData2.rawSignal - repmat(bg2, size(lData2.rawSignal, 1), size(lData2.rawSignal, 2), 1));
mTime2 = lData2.mTime;
if flagOLCor
    ol = load(olFile);
    signal = signal ./ repmat(transpose(ol.ov), 1, size(lData2.rawSignal, 2));
end
rcs2 = signal2 .* repmat(range2, 1, size(lData2.rawSignal, 3)).^2;   % 距离修正信号
snr2 = signal2 ./ repmat(transpose(squeeze(noise2)), size(lData2.rawSignal, 2), 1);   % 信噪比

%% correlation calculation
Rarr = NaN(1, size(signal1, 2));
for iPrf = 1:size(signal1, 2)
    Rarr(iPrf) = corr(signal1(:, iPrf), signal2(:, iPrf), 'type', 'Pearson');
end

%% Display (THI)
figure('Name', 'visibility lidar intercomparison', 'Position', [100, 100, 600, 400], 'Color', 'w');

figpos = subfigPos([0.1, 0.14, 0.78, 0.77], 2, 1, 0, 0.04);

subplot('Position', figpos(1, :), 'Units', 'Normalized');

hold on;
p1 = pcolor(mTime1, range1 * 1e-3, rcs1); hold on;
p1.EdgeColor = 'None';
hold off;

ylabel('Height (km)');
title(sprintf('%s', datestr(mTime1(1), 'yyyy-mm-dd')), 'FontSize', 12, 'FontWeight', 'Bold', 'Interpreter', 'none');

set(gca, 'XTickLabel', '', 'XMinorTick', 'on', 'XTick', mTime1(1):datenum(0, 1, 0, 3, 0, 0):mTime1(end), 'YMinorTick', 'on', ...
        'Box', 'on', 'LineWidth', 1, 'TickDir', 'out', 'FontSize', 12, 'layer', 'top');
ax1 = gca;
ax1.XAxis.MinorTickValues = min(mTime1):datenum(0, 1, 0, 1, 0, 0):max(mTime1);

text(0.02, 0.88, '(a) VIS 1', 'Units', 'Normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'Bold');

xlim([min(mTime1), max(mTime1)]);
ylim(hRange);
caxis([0, 3e10]);
colormap(gca, 'jet');

cb = colorbar('Position', [figpos(1, 1) + figpos(1, 3) + 0.03, figpos(1, 2) + 0.05, 0.02, figpos(1, 4) - 0.10], 'Units', 'Normalized');
set(gca, 'TickDir', 'out', 'Box', 'on');
titleHandle = get(cb, 'Title');
set(titleHandle, 'string', '[a.u.]', 'FontSize', 12);

subplot('Position', figpos(2, :), 'Units', 'Normalized');

hold on;
p1 = pcolor(mTime2, range2 * 1e-3, rcs2); hold on;
p1.EdgeColor = 'None';
hold off;

ylabel('Height (km)');
xlabel('Time (LT)');

set(gca, 'XMinorTick', 'on', 'XTick', mTime2(1):datenum(0, 1, 0, 3, 0, 0):mTime2(end), 'YMinorTick', 'on', ...
        'Box', 'on', 'LineWidth', 1, 'TickDir', 'out', 'FontSize', 12, 'layer', 'top');
ax1 = gca;
ax1.XAxis.MinorTickValues = min(mTime2):datenum(0, 1, 0, 1, 0, 0):max(mTime2);
datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

text(0.02, 0.88, '(b) VIS2', 'Units', 'Normalized', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'Bold');

xlim([min(mTime2), max(mTime2)]);
ylim(hRange);
caxis([0, 8e10]);
colormap(gca, 'jet');

cb = colorbar('Position', [figpos(2, 1) + figpos(2, 3) + 0.03, figpos(2, 2) + 0.05, 0.02, figpos(2, 4) - 0.10], 'Units', 'Normalized');
set(gca, 'TickDir', 'out', 'Box', 'on');
titleHandle = get(cb, 'Title');
set(titleHandle, 'string', '[a.u.]', 'FontSize', 12);

%% correlation plot
figure('Name', 'correlation between two lidars', 'Position', [100, 100, 500, 250], 'Color', 'w');

hold on;
plot(mTime1, Rarr, 'b', 'LineWidth', 2);
yline(0.9, '--r');
yline(1.0, '--k');
hold off;
xlabel('Height (km)');
ylabel('Correlation Coefficient');
xlim([min(mTime2), max(mTime2)]);
ylim([0.5, 1.1]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'LineWidth', 1, 'TickDir', 'out', 'FontSize', 12);
datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');