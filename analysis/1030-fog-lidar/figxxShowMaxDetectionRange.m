clc;
close all;
global LISMO_VARS;

%% 参数设置
dataFolder = 'D:\backup\Data\1030-fog-lidar\data\1030\2023-11-07';
deadtime = 0;
firstRangeBin = 15;
%tRange = [datenum(2023, 8, 2, 20, 0, 0), datenum(2023, 8, 3, 5, 0, 0)];
 tRange = [datenum(2023, 11, 7, 0, 0, 0), datenum(2023, 11, 7, 11, 20, 0)];

%% 读取数据
data = readALADats(dataFolder, 'tRange', tRange);

%% 数据预处理
rawSigPCR = data.rawSignal / (50 * data.nShots(1) * 1e-3);
corSigPCR = rawSigPCR ./ (1 - deadtime * rawSigPCR * 1e-3);
corSigPC = corSigPCR * 50 * data.nShots(1) * 1e-3;
bg = nanmean(corSigPC(:, (end - 80):(end - 10), :), 2);
sigPC = corSigPC - repmat(bg, 1, size(corSigPC, 2), 1);
height = ((1:size(sigPC, 2)) - firstRangeBin + 0.5) * data.hRes(1);
rcs = sigPC .* repmat(reshape(height, 1, length(height), 1), size(sigPC, 1), 1, size(sigPC, 3)).^2;

%% 数据可视化

% 连续观测时空高度图
figure('Position', [0, 0, 600, 250], 'Units', 'Pixels', 'Color', 'w');

p1 = pcolor(data.mTime, height / 1e3, squeeze(rcs(2, :, :) / 1e9));
p1.EdgeColor = 'None';

xlabel('Local Time');
ylabel('Distance (km)');
title(sprintf('Range corrected signal (far-range channel)'));
text(0, -0.2, sprintf('%s', datestr(data.mTime(1), 'yyyy-mm-dd')), 'FontWeight', 'Bold', 'Units', 'Normalized');

caxis([0, 7]);
colormap('jet');

set(gca, 'TickDir', 'out', 'YMinorTick', 'on', 'Box', 'on');
xlim([min(data.mTime), max(data.mTime)]);
ylim([0, 20]);
datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

colorbar();

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'figxx-continuous-measurement.png'), '-r300');

% 连续观测时空高度图
figure('Position', [0, 0, 600, 150], 'Units', 'Pixels', 'Color', 'w');

p1 = pcolor(data.mTime, height, squeeze(rcs(1, :, :) / 1e9));
p1.EdgeColor = 'None';

xlabel('Local Time');
ylabel('Distance (m)');
title(sprintf('Range corrected signal (near-range channel)'));
text(0, -0.4, sprintf('%s', datestr(data.mTime(1), 'yyyy-mm-dd')), 'FontWeight', 'Bold', 'Units', 'Normalized');

caxis([0, 0.1]);
colormap('jet');

set(gca, 'TickDir', 'out', 'YMinorTick', 'on', 'Box', 'on', 'YTick', 0:50:200);
xlim([min(data.mTime), max(data.mTime)]);
ylim([0, 200]);
datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

colorbar();

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'figxx-continuous-measurement-near-range.png'), '-r300');

%% 最远探测距离展示图
figure('Position', [0, 10, 500, 270], 'Units', 'Pixels', 'Color', 'w');

% tForDetectionRange = [datenum(2023, 8, 3, 0, 0, 0), datenum(2023, 8, 3, 0, 30, 0)];
tForDetectionRange = [datenum(2023, 11, 7, 11, 10, 0), datenum(2023, 11, 7, 11, 30, 0)];
isInAVGRange = (data.mTime >= tForDetectionRange(1)) & (data.mTime <= tForDetectionRange(2));
snr = smooth(nansum(sigPC(2, :, isInAVGRange), 3), 4) ./ sqrt(smooth(nansum(corSigPC(2, :, isInAVGRange), 3), 4)) * 2;
startSearchIdx = 200;
idxSNRlt3 = find((snr(startSearchIdx:end) < 3) & (~ isnan(snr(startSearchIdx:end))), 1, 'first') + startSearchIdx;

snr(snr <= 0) = NaN;
p1 = semilogy(height / 1e3, snr, '-b', 'LineWidth', 2);  hold on;
p2 = semilogy([0, 100], [3, 3], '--r');
p3 = scatter(height(idxSNRlt3) / 1e3, snr(idxSNRlt3 - 2), 15, 'Marker', 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

text(height(idxSNRlt3) / 1e3 + 0.5, snr(idxSNRlt3 + 3)*10, sprintf('Distance: %5.2fkm', height(idxSNRlt3) / 1e3), 'Units', 'Data', 'FontSize', 12, 'FontWeight', 'bold');

xlabel('Distance (km)');
ylabel('SNR');
title(sprintf('Maximum detection range (far-range channel)'));

xlim([0, 20]);
ylim([1e-1, 1e4]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.03, 0.02], 'FontSize', 12);

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'figxxMaximumdetectionrange-fr.png'), '-r300');

figure('Position', [0, 10, 500, 270], 'Units', 'Pixels', 'Color', 'w');
isInAVGRange = (data.mTime >= tForDetectionRange(1)) & (data.mTime <= tForDetectionRange(2));
snr = nansum(sigPC(1, :, isInAVGRange), 3) ./ sqrt(nansum(corSigPC(1, :, isInAVGRange), 3));
% startSearchIdx = 200;
% idxSNRlt3 = find((snr(startSearchIdx:end) < 3) & (~ isnan(snr(startSearchIdx:end))), 1, 'first') + startSearchIdx;

snr(snr <= 0) = NaN;
p1 = semilogy(height, snr, '-b', 'LineWidth', 2);  hold on;
% p2 = semilogy([0, 100], [3, 3], '--r');
% p3 = scatter(height(idxSNRlt3) / 1e3, snr(idxSNRlt3 - 2), 15, 'Marker', 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

% text(height(idxSNRlt3) / 1e3 + 0.5, snr(idxSNRlt3 + 3)*10, sprintf('Distance: %5.2fkm', height(idxSNRlt3) / 1e3), 'Units', 'Data', 'FontSize', 12, 'FontWeight', 'bold');

xlabel('Distance (m)');
ylabel('SNR');
title(sprintf('Near-range channel'));

xlim([0, 500]);
ylim([1e-1, 1e4]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.03, 0.02], 'FontSize', 12);

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'figxxMaximumdetectionrange-nr.png'), '-r300');