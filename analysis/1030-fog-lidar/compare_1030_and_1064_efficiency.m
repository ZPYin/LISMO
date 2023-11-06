% 1064 vs 1030 efficiency comparison
% The experiment was done by Qianyuan Chen, at Oct. 11, 2023
% The data was saved at 2023-09-28

clc;
close all;

global LISMO_VARS;

%% Parameter Definition
dataPath = 'C:\Users\ZPYin\Documents\Data\1064nm_cqy\2023-09-28';
tRange1030 = [datenum(2023, 9, 28, 20, 46, 0), datenum(2023, 9, 28, 20, 58, 0)];
tRange1064 = [datenum(2023, 9, 28, 19, 42, 0), datenum(2023, 9, 28, 19, 52, 0)];
deadtime1030 = 22;
deadtime1064 = 22;
firstRangeBin1030 = 15;
firstRangeBin1064 = 15;

data1030 = readALADats(dataPath, 'tRange', tRange1030);
data1064 = readALADats(dataPath, 'tRange', tRange1064);

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

heightRange = [1000, 6000];
isIn1030 = (height1030 >= heightRange(1)) & (height1030 <= heightRange(2));
isIn1064 = (height1064 >= heightRange(1)) & (height1064 <= heightRange(2));
total1030 = sum(sum(rcs1030(2, isIn1030, :), 3), 2);
total1064 = sum(sum(rcs1064(1, isIn1064, :), 3), 2);
efficiencyRatio = total1030 / total1064;

%% Data Display

% profile with full range
figure('Position', [0, 0, 500, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(height1030 / 1e3, squeeze(sum(rcs1030(2, :, :), 3)), '-', 'Color', [199, 101, 74] / 255, 'LineWidth', 2, 'DisplayName', '1030 nm');
p2 = plot(height1064 / 1e3, squeeze(sum(rcs1064(1, :, :), 3)), '-', 'Color', [54, 63, 69] / 255, 'LineWidth', 2, 'DisplayName', '1064 nm');

p3 = patch([height1030 / 1e3, fliplr(height1030 / 1e3)], [squeeze(sum(rcs1030(2, :, :), 3)), fliplr(zeros(1, length(height1030)))], [217, 163, 120] / 255, 'facealpha', 0.5);
p4 = patch([height1064 / 1e3, fliplr(height1064 / 1e3)], [squeeze(sum(rcs1064(1, :, :), 3)), fliplr(zeros(1, length(height1064)))], [94, 124, 136] / 255, 'facealpha', 0.5);

xlim([0, 6]);
ylim([0, 2e11]);

xlabel('距离 (千米)');
ylabel('距离修正信号 (a.u.)');

title(sprintf('实验时间: %s', datestr(data1030.mTime(1), 'yyyy-mm-dd')));

text(0.45, 0.55, sprintf('Received Energy 1030: %4.2e\nReceived Energy 1064: %4.2e\nLaser Power 1030: ~0.6 W\nLaser Power 1064: ~0.3W \n', total1030, total1064), 'Units', 'Normalized', 'FontWeight', 'Bold');

text(0.48, 0.34, sprintf('\\eta_{1030/1064} = %3.1f', total1030 / total1064 * 0.3 / 0.6), 'Units', 'Normalized', 'FontWeight', 'Bold', 'Color', 'r', 'FontSize', 14);

set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'linear', 'TickLength', [0.03, 0.03]);

legend([p1, p2], 'Location', 'NorthEast');

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'efficiency-comparison_1030-1064.png'), '-r300');