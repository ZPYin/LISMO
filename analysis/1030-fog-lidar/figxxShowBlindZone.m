clc;
global LISMO_VARS;

%% 参数设置
dataFolder = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\data\1030\2023-08-04';
deadtime = 22;
firstRangeBin = 16;
tRange = [datenum(2023, 8, 4, 23, 0, 0), datenum(2023, 8, 5, 4, 0, 0)];
visFiles = {'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\visibility-sensor\PreVisi[2023-08-04].csv', ...
            'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\visibility-sensor\PreVisi[2023-08-05].csv'};

%% 读取数据
data = readALADats(dataFolder, 'tRange', tRange);

% 能见度数据
visTime = [];
vis10min = [];
dataVis = readV35(visFiles{1});
visTime = dataVis.mTime;
vis10min = dataVis.vis10min;
dataVis = readV35(visFiles{2});
visTime = cat(2, visTime, dataVis.mTime);
vis10min = cat(2, vis10min, dataVis.vis10min);

%% 数据预处理
rawSigPCR = data.rawSignal / (50 * data.nShots(1) * 1e-3);
corSigPCR = rawSigPCR ./ (1 - deadtime * rawSigPCR * 1e-3);
corSigPC = corSigPCR * 50 * data.nShots(1) * 1e-3;
bg = nanmean(corSigPC(:, (end - 80):(end - 10), :), 2);
sigPC = corSigPC - repmat(bg, 1, size(corSigPC, 2), 1);
height = ((1:size(sigPC, 2)) - firstRangeBin + 0.5) * data.hRes(1);
rcs = sigPC .* repmat(reshape(height, 1, length(height), 1), size(sigPC, 1), 1, size(sigPC, 3)).^2;

%% 数据可视化

%% 连续观测时空高度图(far-range)
figure('Position', [0, 0, 500, 500], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', [0.1, 0.55, 0.82, 0.37], 'Units', 'normalized');
rcsLog = squeeze(rcs(2, :, :) / 1e9);
rcsLog(rcsLog <= 0) = NaN;
p1 = pcolor(data.mTime, height / 1e3, log10(rcsLog));
p1.EdgeColor = 'None';

xlabel('');
ylabel('Distance (km)');
title(sprintf('Logarithm of range corrected signal'), 'FontSize', 11, 'FontWeight', 'bold');
text(0.9, 0.9, '(a)', 'Units', 'Normalized', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'r');

caxis([-0.5, 1]);
colormap('jet');

set(gca, 'TickDir', 'out', 'YMinorTick', 'on', 'Box', 'on', 'XTick', linspace(min(data.mTime), max(data.mTime), 6), 'XTickLabel', '');
xlim([min(data.mTime), max(data.mTime)]);
ylim([0, 20]);

colorbar('Position', [0.93, 0.6, 0.02, 0.27], 'Units', 'normalized');

subplot('Position', [0.1, 0.29, 0.82, 0.23], 'Units', 'Normalized');

rcsLog = squeeze(rcs(1, :, :) / 1e9);
rcsLog(rcsLog <= 0) = NaN;
p1 = pcolor(data.mTime, height, log10(rcsLog));
p1.EdgeColor = 'None';

xlabel('');
ylabel('Distance (m)');

caxis([-2, -1]);
colormap('jet');

set(gca, 'TickDir', 'out', 'YMinorTick', 'on', 'Box', 'on', 'YTick', 0:50:200, 'XTick', linspace(min(data.mTime), max(data.mTime), 6), 'XTickLabel', '');
xlim([min(data.mTime), max(data.mTime)]);
ylim([0, 200]);
text(0.9, 0.9, '(b)', 'Units', 'Normalized', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'r');

colorbar('Position', [0.93, 0.32, 0.02, 0.17], 'Units', 'Normalized');

subplot('Position', [0.1, 0.1, 0.82, 0.15], 'Units', 'Normalized');

plot(visTime, vis10min / 1e3, '-', 'Color', 'b'); hold on;

xlabel('Local Time');
ylabel('Visibility (km)');

set(gca, 'TickDir', 'out', 'YMinorTick', 'on', 'Box', 'on', 'XTick', linspace(min(data.mTime), max(data.mTime), 6));
xlim([min(data.mTime), max(data.mTime)]);
ylim([2e1, 5e1]);
datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');
text(-0.1, -0.5, datestr(data.mTime(1), 'yyyy-mm-dd'), 'Units', 'normalized', 'FontWeight', 'bold');

text(0.9, 0.86, '(c)', 'Units', 'Normalized', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'r');

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'blindzone-demo.png'), '-r300');

%% profiles
figure('Position', [0, 10, 500, 250], 'Units', 'Pixels', 'Color', 'w');

rcsFR = squeeze(nansum(sigPC(2, :, :), 3));
rcsNR = squeeze(nansum(sigPC(1, :, :), 3));

p1 = plot(height /30*200 / 1e6, rawSigPCR(2,:,1), '-b', 'LineWidth', 2, 'DisplayName', 'Far-range');  hold on;
p2 = plot(height/30*200 / 1e6, rawSigPCR(1,:,1), '-g', 'LineWidth', 2, 'DisplayName', 'Near-range');  hold on;

xlabel('time (\mus)');
ylabel('Signal (photon count)');

xlim([0, 30]);
ylim([1, 1e7]);

l1 = legend([p1, p2], 'Location', 'NorthEast');
l1.FontSize = 11;

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.03, 0.02]);

export_fig(gcf, 'signal-blindzone.png', '-r300');

figure('Position', [0, 10, 500, 250], 'Units', 'Pixels', 'Color', 'w');

rcsFR = squeeze(nansum(rcs(2, :, :), 3));
rcsNR = squeeze(nansum(rcs(1, :, :), 3));

p1 = semilogy(height, rcsFR, '-b', 'LineWidth', 2, 'DisplayName', 'Far-range');  hold on;
p2 = semilogy(height, rcsNR * 3e1, '-g', 'LineWidth', 2, 'DisplayName', 'Near-range * 30');  hold on;

xlabel('Distance (m)');
ylabel('Range corrected signal (a.u.)');

xlim([0, 1500]);
ylim([1e6, 1e12]);

l1 = legend([p1, p2], 'Location', 'NorthEast');
l1.FontSize = 11;
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.03, 0.02]);

export_fig(gcf, 'rcs-blindzone.png', '-r300');