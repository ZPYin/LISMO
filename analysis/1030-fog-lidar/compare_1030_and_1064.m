clc;
close all;
global LISMO_VARS;

%% 参数设置
dataFolder1030 = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\data\1030\2023-07-25';
dataFolder1064 = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\1064data\2023-07-25';
deadtime1030 = 22;
deadtime1064 = 22;
firstRangeBin1030 = 15;
firstRangeBin1064 = 15;
tRange = [datenum(2023, 7, 26, 0, 30, 0), datenum(2023, 7, 26, 0, 50, 0)];

%% 读取数据
data1030 = readALADats(dataFolder1030, 'tRange', tRange);
data1064 = readALADats(dataFolder1064, 'tRange', tRange);

%% 数据预处理

% 1030 nm
rawSigPCR1030 = data1030.rawSignal / (50 * data1030.nShots(1) * 1e-3);
corSigPCR1030 = rawSigPCR1030 ./ (1 - deadtime1030 * rawSigPCR1030 * 1e-3);
corSigPC1030 = corSigPCR1030 * 50 * data1030.nShots(1) * 1e-3;
bg1030 = nanmean(corSigPC1030(:, (end - 80):(end - 10), :), 2);
sigPC1030 = corSigPC1030 - repmat(bg1030, 1, size(corSigPC1030, 2), 1);
height1030 = ((1:size(sigPC1030, 2)) - firstRangeBin1030 + 0.5) * data1030.hRes(1);
rcs1030 = sigPC1030 .* repmat(reshape(height1030, 1, length(height1030), 1), size(sigPC1030, 1), 1, size(sigPC1030, 3)).^2;

% 1064 nm
rawSigPCR1064 = data1064.rawSignal / (50 * data1064.nShots(1) * 1e-3);
corSigPCR1064 = rawSigPCR1064 ./ (1 - deadtime1064 * rawSigPCR1064 * 1e-3);
corSigPC1064 = corSigPCR1064 * 50 * data1064.nShots(1) * 1e-3;
bg1064 = nanmean(corSigPC1064(:, (end - 80):(end - 10), :), 2);
sigPC1064 = corSigPC1064 - repmat(bg1064, 1, size(corSigPC1064, 2), 1);
height1064 = ((1:size(sigPC1064, 2)) - firstRangeBin1064 + 0.5) * data1064.hRes(1);
rcs1064 = sigPC1064 .* repmat(reshape(height1064, 1, length(height1064), 1), size(sigPC1064, 1), 1, size(sigPC1064, 3)).^2;

%% Data Display

% profile with full range
figure('Position', [0, 0, 500, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(height1030, squeeze(sum(rcs1030(2, :, :), 3)), '-', 'Color', [15, 15, 172] / 255, 'LineWidth', 2, 'DisplayName', 'Far-Range');
p2 = plot(height1030, squeeze(sum(rcs1030(1, :, :), 3)), '-', 'Color', [56, 87, 35] / 255, 'LineWidth', 2, 'DisplayName', 'Near-Range');
p3 = plot(height1064, squeeze(sum(rcs1064(1, :, :), 3)), '-', 'Color', 'r', 'LineWidth', 1, 'DisplayName', '1064');

xlim([1, 20000]);

xlabel('距离 (米)');
ylabel('距离修正信号');

title(sprintf('累计时间: %s - %s', datestr(data1030.mTime(1), 'yyyy-mm-dd HH:MM'), datestr(data1030.mTime(end), 'HH:MM')));

set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log', 'TickLength', [0.03, 0.03]);

legend([p1, p2, p3], 'Location', 'NorthEast');

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'test-measurement-single-profile_1030_vs_1064.png'), '-r300');

% % profile with limited range
% figure('Position', [0, 0, 500, 300], 'Units', 'Pixels', 'Color', 'w');

% hold on;
% p1 = plot(1:size(corSignal1030, 1), squeeze(corSignal1030(:, 2, iProf)), '-', 'Color', [15, 15, 172] / 255, 'LineWidth', 2, 'DisplayName', 'Far-Range');
% p2 = plot(1:size(corSignal1030, 1), squeeze(corSignal1030(:, 1, iProf)), '-', 'Color', [56, 87, 35] / 255, 'LineWidth', 2, 'DisplayName', 'Near-Range');

% xlim([1, 200]);

% xlabel('Range Bin');
% ylabel('Raw Signal (MHz)');

% title(sprintf('Time: %s', datestr(data.mTime(iProf), 'yyyy-mm-dd HH:MM:SS')));

% set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log');

% %legend([p1, p2], 'Location', 'NorthEast');

% export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'test-measurement-single-profile-zoomin.png'), '-r300');

%% snr
figure('Position', [0, 0, 500, 300], 'Units', 'Pixels', 'Color', 'w');
hold on;
snrFR = sum(sigPC1030(2, :, :), 3) ./ sqrt(sum(corSigPC1030(2, :, :), 3));
snrNR = sum(sigPC1030(1, :, :), 3) ./ sqrt(sum(corSigPC1030(1, :, :), 3));
snr1064 = sum(sigPC1064(1, :, :), 3) ./ sqrt(sum(corSigPC1064(1, :, :), 3));
p1 = semilogy(height1030, snrFR, '-', 'Color', [15, 15, 172] / 255, 'LineWidth', 2, 'DisplayName', 'Far-Range');
p2 = semilogy(height1030, snrNR, '-', 'Color', [56, 87, 35] / 255, 'LineWidth', 2, 'DisplayName', 'Near-Range');
p3 = semilogy(height1064, snr1064, '-', 'Color', 'r', 'LineWidth', 1, 'DisplayName', '1064');

p4 = semilogy([0, 1e5], [3, 3], '--r');
xlim([0, 20000]);

xlabel('探测距离 (米)');
ylabel('信噪比');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log', 'TickLength', [0.03, 0.03]);

legend([p1, p2, p3], 'Location', 'NorthEast');

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'test-measurement-SNR_1030_vs_1064.png'), '-r300');

% %% range corrected signal
% figure('Position', [0, 0, 500, 300], 'Units', 'Pixels', 'Color', 'w');

% hold on;
% p1 = plot(height, smooth(nanmean(squeeze(rcs(:, 2, iProf)), 2), 4), '-', 'Color', [15, 15, 172] / 255, 'LineWidth', 2, 'DisplayName', 'Far-Range');
% p2 = plot(height, nanmean(squeeze(rcs(:, 1, iProf)), 2), '-', 'Color', [56, 87, 35] / 255, 'LineWidth', 2, 'DisplayName', 'Near-Range');

% xlim([0, max(height)]);

% xlabel('Distance [m]');
% ylabel('range cor. sig. [a.u.]');

% set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'linear', 'TickLength', [0.03, 0.03]);

% legend([p1, p2], 'Location', 'NorthEast');

% export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'test-measurement-single-profile-rcs.png'), '-r300');

% %% range corrected signal (near-range)
% figure('Position', [0, 0, 500, 300], 'Units', 'Pixels', 'Color', 'w');

% hold on;
% p1 = plot(height, nanmean(squeeze(rcs(:, 2, (iProf - 60):iProf)), 2), '-', 'Color', [15, 15, 172] / 255, 'LineWidth', 2, 'DisplayName', 'Far-Range');
% p2 = plot(height, nanmean(squeeze(rcs(:, 1, (iProf - 60):iProf)), 2), '--', 'Color', [56, 87, 35] / 255, 'LineWidth', 1, 'DisplayName', 'Near-Range');
% p3 = plot(height, nanmean(squeeze(rcs(:, 1, (iProf - 60):iProf)), 2) * 140000.0, '-', 'Color', [56, 87, 35] / 255, 'Marker', 'o', 'MarkerSize', 5, 'LineWidth', 2, 'DisplayName', 'Near-Range (x k)');

% xlim([0, 500]);

% xlabel('Distance [m]');
% ylabel('range cor. sig. [a.u.]');

% title(sprintf('Time: %s', datestr(data.mTime(iProf), 'yyyy-mm-dd HH:MM:SS')));

% set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log', 'TickLength', [0.03, 0.03]);

% legend([p1, p2, p3], 'Location', 'SouthEast');

% export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'test-measurement-single-profile-rcs-near-range.png'), '-r300');