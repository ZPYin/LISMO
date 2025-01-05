clc;
% close all;

%% Parameter Definition
tRange = datenum(2024, 11, 1, 0, 20, 0);   % 选择的廓线时间
visSensorFile = 'vis-sensor-data.mat';   % 前向散射能见度仪数据文件
l0Folder = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\lidar';   % 雷达原始数据文件目录
l1Folder = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\results';   % 雷达产品文件目录
saveFolder = 'C:\Users\ZPYin\Desktop';   % 输出结果目录
overlapFile = 'lk_20241101_overlap_factor.txt';   % 重叠因子文件，重叠因子通过水平方法计算
flagOLCor = false;   % 是否进行重叠因子修正
lr = 50;   % 雷达比
dist2_998 = 4.5;   % I0998站点与激光雷达距离（千米）
dist2_297 = 3.2;   % I0297站点与激光雷达距离（千米）
refH = [8, 9];   % Fernald方法参考高度（千米）

%% Single Profile
filename = sprintf('%s.pb', datestr(tRange, 'yyyymmddHHMMSS'));

% read lidar data
thisDate = datenum(filename(1:(end - 3)), 'yyyymmddHHMMSS');
thisData = readVisLidarL0(fullfile(l0Folder, datestr(thisDate, 'yyyy-mm-dd'), filename));

% read level 1
thisL1 = readVisLidarL1(fullfile(l1Folder, datestr(thisDate, 'yyyy-mm-dd'), filename));

% read vis sensor
load(visSensorFile);

% read overlap factor
ol = [];
if exist(overlapFile, 'file') == 2
    fid = fopen(overlapFile, 'r');

    tmp = textscan(fid, '%f%f', 'HeaderLines', 1, 'Delimiter', ' ');
    olHeight = tmp{1};
    olVal = tmp{2};

    fclose(fid);
end

% find pretrigger
% figure;
% plot(thisData.sig);

% xlabel('Range Bin');
% ylabel('Signal (a.u.)');

% xlim([0, 200]);

% set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

% preprocessing
nPretrigger = 12;
thisBG = nanmean(thisData.sig((end - 200):end));
thisSigCor = thisData.sig(nPretrigger:end) - thisBG;
height = transpose(((1:length(thisSigCor)) - 0.5) * thisData.hRes * 1e-3);
ol = interp1(olHeight, olVal, height * 1e3);
if flagOLCor
    thisSigCor = thisSigCor ./ ol;
end
thisRCS = thisSigCor .* height.^2;

% Rayleigh scattering
[mBsc, mExt] = MolModel(ones(size(height)) * 0.005 * 1e3, 1064, 'meteor', 'standard_atmosphere');
[~, mExt550] = MolModel(ones(size(height)) * 0.005 * 1e3, 550, 'meteor', 'standard_atmosphere');
mAttn = mBsc .* exp(-2 * nancumsum(mExt .* [height(1); diff(height)] * 1e3));
ratioL2M = nansum(thisRCS(1000:1200)) / nansum(mAttn(1000:1200));

% display
figure('Position', [0, 20, 600, 250], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(height, thisRCS ./ ratioL2M, '-k', 'DisplayName', '激光雷达');
p2 = plot(height, mAttn, '--r', 'DisplayName', '大气分子');
hold off;

xlabel('距离 (千米)');
ylabel('距离修正信号');

xlim([0, 30]);
ylim([1e-8, 1e-3]);

set(gca, 'YScale', 'log', 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

legend([p1, p2], 'Location', 'NorthEast');

% export_fig(gcf, fullfile(saveFolder, sprintf('%s_rcs.png', datestr(tRange(1), 'yyyy-mm-dd'))), '-r300');

% lidar calibration
isInHCaliRange = (height >= refH(1)) & (height <= refH(2));
try
    [~, extRef] = chi2fit(height(isInHCaliRange) * 1e3, log(thisRCS(isInHCaliRange)), 0.1 * ones(size(thisRCS(isInHCaliRange))));
catch
    extRef = 1e-4;
end
extRef = extRef / -2;
aBsc = transpose(fernald(height * 1e3, thisSigCor, thisBG, lr, refH * 1e3, extRef / lr, mBsc, 4));
aBscCor = aBsc;
aBscCor(height <= 1.5) = aBscCor(find(height >= 1.5, 1));
lc = thisRCS ./ ((aBscCor + mBsc) .* exp(-2 * nancumsum(aBscCor * lr + mExt).* [height(1); diff(height)] * 1e3));
lcMean = nanmean(lc(500:600));

% extinction retrieval
[~, ext_Holger] = extRet_Holger(height' * 1e3, thisSigCor', ...
                           'calibration_constant', lcMean * 1e6, ...
                           'fullOverlapR', 1200, ...
                           'elevation_angle', 0);
ext_Holger = transpose(ext_Holger);
ext_Xian = extRet_Xian(height * 1e3, thisSigCor * 3000 * 50 * 1e-3, thisBG * 3000 * 50 * 1e-3, 'minSNR', 0.2, 'rangeFullOverlap', 1200);
ext_Fernald = extRet_Fernald(height * 1e3, thisSigCor * 3000 * 50 * 1e-3, thisBG * 3000 * 50 * 1e-3, mBsc, mExt, 'minSNR', 3, 'hFullOL', 1200);
ext_LK_1 = interp1(thisL1.height, thisL1.extinction, height * 1e3);

% visiblity
AEConvFactor = 2.52;   % 波长指数为1.4
vis_Holger = ext2vis(ext_Holger * AEConvFactor + mExt550);
vis_Xian = ext2vis(ext_Xian);
vis_Fernald = ext2vis(ext_Fernald * AEConvFactor + mExt550);
vis_LK_1 = interp1(thisL1.height, thisL1.vis, height * 1e3);

% find vis sensor data
refIdx998 = find(height >= dist2_998, 1);
[minDiff, visIdx998] = min(abs(vis998.mTime - thisDate));
thisVisRef998 = NaN;
if minDiff <= datenum(0, 1, 0, 0, 5, 0)
    thisVisRef998 = vis998.vis(visIdx998);
end
vis998_Holger = vis_Holger(refIdx998);
vis998_Xian = vis_Xian(refIdx998);
vis998_Fernald = vis_Fernald(refIdx998);
vis998_LK = vis_LK_1(refIdx998);

refIdx297 = find(height >= dist2_297, 1);
[minDiff, visIdx297] = min(abs(vis297.mTime - thisDate));
thisVisRef297 = NaN;
if minDiff <= datenum(0, 1, 0, 0, 5, 0)
    thisVisRef297 = vis297.vis(visIdx297);
end
vis297_Holger = vis_Holger(refIdx297);
vis297_Xian = vis_Xian(refIdx297);
vis297_Fernald = vis_Fernald(refIdx297);
vis297_LK = vis_LK_1(refIdx297);

%% display profile

% extinction
figure('Position', [0, 20, 600, 250], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(height, ext_Fernald * 1e6, '-k', 'DisplayName', 'Fernald算法');
plot(height(isInHCaliRange), ext_Fernald(isInHCaliRange) * 1e6, '-r', 'LineWidth', 2);
p2 = plot(height, ext_LK_1 * 1e3, '-r', 'DisplayName', '蓝科光电');
p3 = plot(height, ext_Holger * 1e6, '-g', 'DisplayName', '前向迭代');
p4 = plot(height, ext_Xian * 1e6 - mExt * 1e6, '-b', 'DisplayName', '大舜算法');
plot([height(refIdx998), height(refIdx998)], [-1e5, 1e5], '--', 'color', 'k');
plot([height(refIdx297), height(refIdx297)], [-1e5, 1e5], '--', 'color', 'k');
hold off;

text(height(refIdx998) * 1.1, ext_Fernald(refIdx998) * 1e6 * 1.1, sprintf('I0998'), 'Units', 'Data', 'FontSize', 11, 'color', 'k');
text(height(refIdx297) * 1.1, ext_Fernald(refIdx297) * 1e6 * 1.1, sprintf('I0297'), 'Units', 'Data', 'FontSize', 11, 'color', 'k');
text(0.3, 0.8, sprintf('参考距离: %5.2f-%5.2f 千米', refH(1), refH(2)), 'Units', 'Normalized', 'FontSize', 11);

xlabel('距离 (千米)');
ylabel('气溶胶消光系数 (Mm-1)');

xlim([0, 20]);
ylim([0, 500]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

legend([p1, p2, p3, p4], 'Location', 'NorthEast');

% export_fig(gcf, fullfile(saveFolder, sprintf('%s_single_profile_extinction.png', datestr(tRange, 'yyyymmddHHMMSS'))), '-r300');

% visibility
figure('Position', [0, 20, 600, 250], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(height, vis_Fernald * 1e-3, '-k', 'DisplayName', 'Fernald算法');
p2 = plot(height, vis_LK_1, '-r', 'DisplayName', '蓝科光电');
p3 = plot(height, vis_Holger * 1e-3, '-g', 'DisplayName', '前向迭代');
p4 = plot(height, vis_Xian * 1e-3, '-b', 'DisplayName', '大舜算法');
p5 = scatter(height(refIdx998), thisVisRef998 * 1e-3, 30, 'Marker', 's', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'DisplayName', 'I0998');
p6 = scatter(height(refIdx297), thisVisRef297 * 1e-3, 30, 'Marker', 's', 'MarkerEdgeColor', 'cyan', 'MarkerFaceColor', 'cyan', 'DisplayName', 'I0297');
hold off;

xlabel('距离 (千米)');
ylabel('能见度 (千米)');

xlim([0, 20]);
ylim([0, 50]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

legend([p1, p2, p3, p4, p5, p6], 'Location', 'NorthEast');

% export_fig(gcf, fullfile(saveFolder, sprintf('%s_single_profile_vis.png', datestr(tRange, 'yyyymmddHHMMSS'))), '-r300');