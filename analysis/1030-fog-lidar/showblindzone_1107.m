clc;
close all;
global LISMO_VARS;

%% 参数设置
dataFolder = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\data\1030\2023-11-07';
deadtime = 22;
firstRangeBin = 15;
% tRange = [datenum(2023, 7, 17, 0, 0, 0), datenum(2023, 7, 17, 0, 30, 0)];
tRange = [datenum(2023, 11, 7, 13, 18, 0), datenum(2023, 11, 7, 13, 18, 50)];

%% 读取数据
data = readALADats(dataFolder, 'tRange', tRange, 'nMaxBin', 2100);

%% 数据预处理
rawSigPCR = data.rawSignal / (50 * data.nShots(1) * 1e-3);
corSigPCR = rawSigPCR ./ (1 - deadtime * rawSigPCR * 1e-3);
corSigPC = corSigPCR * 50 * data.nShots(1) * 1e-3;
bg = nanmean(corSigPC(:, (end - 80):(end - 10), :), 2);
sigPC = corSigPC - repmat(bg, 1, size(corSigPC, 2), 1);
height = ((1:size(sigPC, 2)) - firstRangeBin + 0.5) * data.hRes(1);
rcs = sigPC .* repmat(reshape(height, 1, length(height), 1), size(sigPC, 1), 1, size(sigPC, 3)).^2;

%% 数据可视化

% snr
figure('Position', [0, 10, 500, 270], 'Units', 'Pixels', 'Color', 'w');

rcsFR = squeeze(nanmean(sigPC(2, :, :), 3));
rcsNR = squeeze(nanmean(sigPC(1, :, :), 3));

p2 = semilogy(height, rcsNR / (50 * data.nShots(1) * 1e-3), '-g', 'LineWidth', 2, 'DisplayName', '近场通道', 'marker', '.', 'markersize', 15);  hold on;
p1 = semilogy(height, rcsFR / (50 * data.nShots(1) * 1e-3), '-b', 'LineWidth', 2, 'DisplayName', '远场通道', 'marker', '.', 'markersize', 15);  hold on;

xlabel('距离 (米)');
ylabel('探测信号 (MHz)');

xlim([0, 800]);
ylim([1e-2, 30]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.03, 0.02], 'layer', 'top');

legend([p1, p2], 'location', 'southeast');

figure('Position', [0, 10, 500, 270], 'Units', 'Pixels', 'Color', 'w');

rcsFR = squeeze(nansum(rcs(2, :, :), 3));
rcsNR = squeeze(nansum(rcs(1, :, :), 3));
ratio = 1;

p1 = semilogy(height, rcsFR, '-b', 'LineWidth', 2, 'DisplayName', '远场通道');  hold on;
p2 = semilogy(height, rcsNR * ratio, '-g', 'LineWidth', 2, 'DisplayName', sprintf('近场通道 * %d', ratio));  hold on;

xlabel('距离 (米)');
ylabel('距离修正信号');

xlim([0, 8000]);
ylim([1e8, 1e11]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.03, 0.02], 'layer', 'top');

legend([p1, p2], 'location', 'northeast');

% snr
figure('Position', [0, 10, 500, 270], 'Units', 'Pixels', 'Color', 'w');

snrFR = squeeze(sigPC(2, :, :) ./ sqrt(corSigPC(2, :, :)));
snrNR = squeeze(sigPC(1, :, :) ./ sqrt(corSigPC(1, :, :)));

heightInterp = 0:7.5:31200;

p2 = semilogy(heightInterp, interp1(height, snrNR, heightInterp), '-g', 'LineWidth', 2, 'DisplayName', '近场通道', 'marker', '.', 'markersize', 10);  hold on;
p1 = semilogy(heightInterp, interp1(height, snrFR, heightInterp), '-b', 'LineWidth', 2, 'DisplayName', '远场通道', 'marker', '.', 'markersize', 10);  hold on;
p3 = plot([0, 1e10], [100, 100], '--r');

xlabel('距离 (米)');
ylabel('信噪比');

xlim([0, 500]);
ylim([1, 1e4]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.03, 0.02], 'layer', 'top');

legend([p1, p2], 'location', 'southeast');

figure('Position', [0, 10, 500, 270], 'Units', 'Pixels', 'Color', 'w');

rcsFR = squeeze(nansum(rcs(2, :, :), 3));
rcsNR = squeeze(nansum(rcs(1, :, :), 3));
ratio = 1;

p1 = semilogy(height, rcsFR, '-b', 'LineWidth', 2, 'DisplayName', '远场通道');  hold on;
p2 = semilogy(height, rcsNR * ratio, '-g', 'LineWidth', 2, 'DisplayName', sprintf('近场通道 * %d', ratio));  hold on;

xlabel('距离 (米)');
ylabel('距离修正信号');

xlim([0, 2000]);
ylim([1e4, 1e11]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.03, 0.02], 'layer', 'top');

legend([p1, p2], 'location', 'northeast');