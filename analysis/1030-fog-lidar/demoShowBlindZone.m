clc;
global LISMO_VARS;

%% 参数设置
dataFolder = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\data\1030\2023-11-07';
deadtime = 22;
firstRangeBin = 16;
% tRange = [datenum(2023, 7, 17, 0, 0, 0), datenum(2023, 7, 17, 0, 30, 0)];
tRange = [datenum(2023, 11, 7, 13, 18, 0), datenum(2023, 11, 7, 13, 19, 0)];

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

% % 连续观测时空高度图
% figure('Position', [0, 0, 450, 250], 'Units', 'Pixels', 'Color', 'w');
% 
% subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
% mTimeInterp = data.mTime(1):datenum(0, 1, 0, 0, 1, 0):data.mTime(end);
% rcsInterp = NaN(length(height), length(mTimeInterp));
% for iTime = 1:length(data.mTime)
%     tIdx = round((data.mTime(iTime) - data.mTime(1)) / datenum(0, 1, 0, 0, 1, 0) + 1);
%     rcsInterp(:, tIdx) = rcs(2, :, iTime);
% end
% p1 = pcolor(mTimeInterp, height / 1e3, rcsInterp / 1e9);
% p1.EdgeColor = 'None';
% 
% xlabel('时间 (小时:分钟)');
% ylabel('距离 (千米)');
% title(sprintf('激光测雾雷达距离修正信号 %s', datestr(data.mTime(1), 'yyyy-mm-dd')));
% 
% caxis([0, 50]);
% colormap('jet');
% 
% set(gca, 'TickDir', 'out', 'YMinorTick', 'on', 'Box', 'on', 'XTick', linspace(data.mTime(1), data.mTime(end), 5), 'layer', 'top');
% xlim([min(data.mTime), max(data.mTime)]);
% ylim([0, 30]);
% datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');
% 
% colorbar('position', [0.9, 0.2, 0.03, 0.6], 'Units', 'Normalized');

% export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', '最小盲区计算使用信号区间展示.png'), '-r300');

%% 最远探测距离展示图
figure('Position', [0, 10, 500, 270], 'Units', 'Pixels', 'Color', 'w');

rcsFR = squeeze(nansum(sigPC(2, :, :), 3));
rcsNR = squeeze(nansum(sigPC(1, :, :), 3));

p2 = semilogy(height / 1e3, rcsNR, '-g', 'LineWidth', 2, 'DisplayName', '近场通道');  hold on;
p1 = semilogy(height / 1e3, rcsFR, '-b', 'LineWidth', 2, 'DisplayName', '远场通道');  hold on;

xlabel('距离 (千米)');
ylabel('探测信号 (光子数)');

xlim([0, 30]);
ylim([1, 1e7]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.03, 0.02], 'layer', 'top');

legend([p1, p2], 'location', 'northeast');

figure('Position', [0, 10, 500, 270], 'Units', 'Pixels', 'Color', 'w');

rcsFR = squeeze(nansum(rcs(2, :, :), 3));
rcsNR = squeeze(nansum(rcs(1, :, :), 3));

p1 = semilogy(height, rcsFR, '-b', 'LineWidth', 2, 'DisplayName', '远场通道');  hold on;
p2 = semilogy(height, rcsNR, '-g', 'LineWidth', 2, 'DisplayName', '近场通道');  hold on;

xlabel('距离 (米)');
ylabel('距离修正信号');

xlim([0, 8000]);
ylim([1e5, 1e12]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.03, 0.02], 'layer', 'top');

legend([p1, p2], 'location', 'northeast');
