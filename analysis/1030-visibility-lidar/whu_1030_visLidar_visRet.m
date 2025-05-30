% 反演得到武汉大学能见度激光雷达产品
% - 气溶胶后向散射系数
% - 气溶胶消光系数
% - 能见度
% 作者：殷振平
% 邮箱：zp.yin@whu.edu.cn
% 日期：2025年5月21日

clc;
close all;

%% Parameter Definition
dataFile = 'C:\Users\zhenp\OneDrive\Desktop\VisibilityInversion (3)\VisibilityInversion (3)\VisibilityInversion\S001-NAP001-Test CMA-017-210927-110120.dat';   % 信号廓线
resPath = 'C:\Users\zhenp\OneDrive\Desktop';   % 产品输出目录
figPath = 'C:\Users\zhenp\OneDrive\Desktop';   % 图片结果输出目录
flagOLCor = false;   % 是否进行重叠因子修正
olFile = 'overlap_20250101.mat';   % 重叠因子文件
AEConvFactor = (1030/550) ^ 1;   % 气溶胶消光系数波长转换因子
hFullOL = 500;   % 完全进入视场高度（m）
lidarRatio = 50;   % 激光雷达比 (sr)
distOffset = -55;   % 预触发点个数（可通过信号第一个峰值进行判断）
visible = 'on';   % 是否进行结果可视化

%% Read Data
if exist(dataFile, 'file')
    thisData = readALADat(dataFile);
else
    warning('Data file does not exist: %s', dataFile);
    return;
end

%% Preprocessing
range = transpose(((1:thisData.nBins) - 0.5 + distOffset) * thisData.hRes(1));
bg = mean(thisData.rawSignal((end - 30):(end - 5), :), 1);
noise = std(thisData.rawSignal((end - 30):(end - 5), :), 0, 1);
signal = thisData.rawSignal - repmat(bg, thisData.nBins, 1);
angle = thisData.ele;
mTime = thisData.mTime;
if flagOLCor
    ol = load(olFile);
    signal = signal ./ repmat(transpose(ol.ov), 1, size(thisData.rawSignal, 2));
end
rcs = signal .* repmat(range, 1, size(thisData.rawSignal, 2)).^2;   % 距离修正信号
snr = signal ./ repmat(noise, thisData.nBins, 1);   % 信噪比
height = range .* sin(angle / 180 * pi);   % 高度（m）

%% Rayleigh Scattering
[mBsc, mExt] = MolModel(height, 1030, 'meteor', 'standard_atmosphere');
[~, mExt550] = MolModel(height, 550, 'meteor', 'standard_atmosphere');

%% extinction retrieval
extMat_Fernald = extRet_Fernald(range, signal, bg, mBsc, mExt, ...
    'snr', snr, ...
    'minSNR', 3, ...
    'hFullOL', 500, ...
    'lr', lidarRatio);
bscMat_Fernald = extMat_Fernald ./ lidarRatio;

%% Extinction to Visibility
visMat_Fernald = ext2vis(extMat_Fernald * AEConvFactor + mExt550);

%% Save Results
fid = fopen(fullfile(resPath, 'test_output.txt'), 'w');

fprintf(fid, '输入信息: \n');
fprintf(fid, '文件名: %s\n', basename(dataFile));

fprintf(fid, '\n反演结果\n');
fprintf(fid, '距离(m); 距离修正信号; 信噪比; 气溶胶后向散射系数 (km-1sr-1); 气溶胶消光系数 (km-1); 能见度 (km);\n');

for iBin = 1:thisData.nBins
    fprintf(fid, '%f; %f; %f; %f; %f; %f;\n', range(iBin), rcs(iBin), snr(iBin), bscMat_Fernald(iBin) * 1e3, extMat_Fernald(iBin) * 1e3, visMat_Fernald(iBin) * 1e-3);
end

fclose(fid);

%% Display

% range corrected signal
figure('Position', [100, 100, 600, 600], 'color', 'w', 'visible', visible);

% range corrected signal
subplot(311);

hold on;
YRcs = rcs;
YRcs(YRcs <= 0) = NaN;
plot(range * 1e-3, YRcs, 'b');
hold off;

xlabel('距离 (km)');
ylabel('距离修正信号 (a.u.)');

xlim([0, 15]);
ylim('auto');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'YTick', 10.^(6:14), ...
    'Box', 'on', 'FontSize', 12, 'YScale', 'log');

% aerosol backscatter and extinction
subplot(312);

hold on;
p1 = plot(range * 1e-3, extMat_Fernald * 1e3, 'b', 'DisplayName', '气溶胶');
p2 = plot(range * 1e-3, mExt * 1e3, 'r', 'DisplayName', '大气分子');
hold off;

xlabel('距离 (km)');
ylabel('消光系数 (km-1)');

xlim([0, 15]);
ylim([-0.005, 0.1]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'YTick', 0:0.2:1, ...
    'Box', 'on', 'FontSize', 12, 'YScale', 'Linear');

l = legend([p1, p2], 'Location', 'NorthEast');
l.FontSize = 11;

% visibility
subplot(313);

hold on;
plot(range * 1e-3, visMat_Fernald * 1e-3, 'b', 'DisplayName', '气溶胶');
plot([-100, 100], [50, 50], '-.k');
hold off;

xlabel('距离 (km)');
ylabel('能见度 (km)');

xlim([0, 15]);
ylim([0, 60]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'YTick', 0:10:50, ...
    'Box', 'on', 'FontSize', 12, 'YScale', 'Linear');

export_fig(gcf, fullfile(figPath, 'test_results_output.png'), '-r300');
