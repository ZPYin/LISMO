%% Parameter Definition
visSensorFile = 'vis-sensor-data.mat';   % 前向散射能见度仪数据文件
l0Folder = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\lidar';   % 雷达原始数据文件目录
l1Folder = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\results';   % 雷达产品文件目录
saveFolder = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\quicklooks';   % 输出结果目录
overlapFile = 'lk_20241101_overlap_factor.txt';   % 重叠因子文件，重叠因子通过水平方法计算
flagOLCor = true;   % 是否进行重叠因子修正
lr = 50;   % 雷达比
dist2_998 = 4.5;   % I0998站点与激光雷达距离（千米）
dist2_297 = 3.2;   % I0297站点与激光雷达距离（千米）
refH = [8, 9];   % Fernald方法参考高度（千米）
nPretrigger = 12;
flagReadData = true;
AEConvFactor = 2.52;   % 波长指数为1.4
visDiffThresh = 2e4;   % [m]

%% Load Data

% vis sensor
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

%% Statistics File
sFid = fopen(fullfile(saveFolder, 'vis_statistics.txt'), 'w');

%% find dates
dateFolders = listdir(l0Folder, '\w*', 1);

for iFolder = 1:1%length(dateFolders)

    close all;

    thisDate = datenum(basename(dateFolders{iFolder}), 'yyyy-mm-dd');
    thisFiles = listfile(dateFolders{iFolder}, '\w*pb', 1);

    filename = basename(thisFiles{1});

    % read lidar data
    thisData = readVisLidarL0(fullfile(dateFolders{iFolder}, filename));

    % read level 1
    thisL1 = readVisLidarL1(fullfile(l1Folder, datestr(thisDate, 'yyyy-mm-dd'), filename));

    % preprocessing
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
    ratioL2M = nansum(thisRCS(2000:2200)) / nansum(mAttn(2000:2200));

    % lidar calibration
    isInHCaliRange = (height >= refH(1)) & (height <= refH(2));
    try
        [~, extRef] = chi2fit(height(isInHCaliRange) * 1e3, log(thisRCS(isInHCaliRange)), 0.1 * ones(size(thisRCS(isInHCaliRange))));
    catch
        extRef = 1e-5;
    end
    extRef = extRef / -2;
    aBsc = transpose(fernald(height * 1e3, thisSigCor, thisBG, lr, refH * 1e3, extRef / lr, mBsc, 4));

    aBscCor = aBsc;
    aBscCor(height <= 1.5) = aBscCor(find(height >= 1.5, 1));
    lc = thisRCS ./ ((aBscCor + mBsc) .* exp(-2 * nancumsum(aBscCor * lr + mExt).* [height(1); diff(height)] * 1e3));
    lcMean = nanmean(lc(500:600)) * 1.3;

    if flagReadData
        % Read Lidar Data
        lidarSig = [];
        laserTemp = [];
        mTime = [];
        ext_LK = [];
        vis_LK = [];
        height_LK = [];
        mTime_LK = [];
        for iFile = 1:length(thisFiles)
            fprintf('Finished %6.2f%%: reading %s\n', (iFile - 1) / length(thisFiles) * 100, thisFiles{iFile});
            thisData = readVisLidarL0(thisFiles{iFile});
            thisData.sig = thisData.sig(nPretrigger:end);

            lidarSig = cat(2, lidarSig, thisData.sig);
            mTime = cat(2, mTime, thisData.startTime);
            laserTemp = cat(2, laserTemp, thisData.laserTemp);

            thisL1 = readVisLidarL1(fullfile(l1Folder, datestr(thisDate, 'yyyy-mm-dd'), basename(thisFiles{iFile})));
            vis_LK = cat(2, vis_LK, thisL1.vis);
            ext_LK = cat(2, ext_LK, thisL1.extinction);
            mTime_LK = cat(2, mTime_LK, thisL1.mTime);
            height_LK = thisL1.height;
        end
    end

    %% signal preprocessing
    thisBG = nanmean(lidarSig((end - 200):end, :), 1);
    lidarSigCor = lidarSig - repmat(thisBG, size(thisBG, 1), 1);
    if flagOLCor
        lidarSigCor = lidarSigCor ./ repmat(ol, 1, size(lidarSigCor, 2));
    end
    lidarRCS = lidarSigCor .* repmat(height, 1, size(lidarSigCor, 2)).^2;

    %% Display

    % rcs
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w');

    hold on;
    p1 = pcolor(mTime, height, lidarRCS);
    p1.EdgeColor = 'none';
    hold off;

    xlabel('时间');
    ylabel('距离 (千米)');
    title('距离修正信号');

    xlim([min(mTime), max(mTime)]);
    ylim([0, 10]);
    caxis([0, 30]);

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'layer', 'Top', 'TickDir', 'out', 'FontSize', 11);

    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    colorbar();
    export_fig(gcf, fullfile(saveFolder, sprintf('%s_RCS.png', datestr(mTime(1), 'yyyy-mm-dd'))), '-r300');

    %% visiblity
    extMat_Xian = NaN(size(lidarSigCor));
    extMat_Holger = NaN(size(lidarSigCor));
    extMat_Fernald = NaN(size(lidarSigCor));
    [MTIMELK, HEIGHTLK] = meshgrid(mTime_LK, height_LK * 1e-3);
    [MTIME, HEIGHT] = meshgrid(mTime, height);
    extMat_LK = interp2(mTime_LK, height_LK, ext_LK, MTIME, HEIGHT);
    % for iPrf = 590:840%size(lidarSigCor, 2)
    for iPrf = 1:size(lidarSigCor, 2)
        thisRCS = lidarRCS(:, iPrf);
        thisRCS(thisRCS <= 0) = NaN;
        isInHCaliRange = (height >= refH(1)) & (height <= refH(2));
        [~, extRef] = chi2fit(height(isInHCaliRange) * 1e3, log(thisRCS(isInHCaliRange)), 0.1 * ones(size(thisRCS(isInHCaliRange))));
        extRef = extRef / -2;

        % extinction retrieval
        ext1 = extRet_Fernald(height * 1e3, lidarSigCor(:, iPrf) * 3000 * 50 * 1e-3, thisBG(iPrf) * 3000 * 50 * 1e-3, mBsc, mExt, 'minSNR', 3, 'hFullOL', 1200, 'calibration_constant', lcMean * 1e6);
        [~, ext2] = extRet_Holger(height' * 1e3, lidarSigCor(:, iPrf)', ...
                                'calibration_constant', lcMean * 1e6, ...
                                'fullOverlapR', 1200, ...
                                'elevation_angle', 0);
        ext3 = extRet_Xian(height * 1e3, lidarSigCor(:, iPrf) * 3000 * 50 * 1e-3, thisBG(iPrf) * 3000 * 50 * 1e-3, 'minSNR', 2, 'rangeFullOverlap', 1200);

        extMat_Xian(:, iPrf) = ext3;
        extMat_Holger(:, iPrf) = ext2;
        extMat_Fernald(:, iPrf) = ext1;
    end

    %%
    visMat_Xian = ext2vis(extMat_Xian);
    visMat_Holger = ext2vis(extMat_Holger * AEConvFactor + mExt550);
    visMat_Fernald = ext2vis(extMat_Fernald * AEConvFactor + mExt550);
    visMat_LK = interp2(MTIMELK, HEIGHTLK, vis_LK * 1e3, MTIME, HEIGHT);

    %% Display

    % Xian
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w');

    subplot('Position', [0.1, 0.2, 0.72, 0.7], 'Units', 'Normalized');

    hold on;
    p1 = pcolor(mTime, height, visMat_Xian * 1e-3);
    p1.EdgeColor = 'none';

    xlabel('时间');
    ylabel('距离 (千米)');
    title(sprintf('%s 能见度 (大舜算法)', datestr(mTime(1), 'yyyy-mm-dd')));

    xlim([min(mTime), max(mTime)]);
    ylim([0, 10]);
    caxis([0, 50]);

    colormap('jet');

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'XTick', linspace(min(mTime), max(mTime), 5), 'layer', 'Top', 'TickDir', 'out', 'FontSize', 11);

    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    cb = colorbar('Position', [0.85, 0.3, 0.03, 0.5], 'Units', 'Normalized');
    titleHandle = get(cb, 'Title');
    set(titleHandle, 'String', '[千米]');
    set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

    export_fig(gcf, fullfile(saveFolder, sprintf('%s_vis_thi_Xian.png', datestr(mTime(1), 'yyyy-mm-dd'))), '-r300');

    %% Holger
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w');

    subplot('Position', [0.1, 0.2, 0.72, 0.7], 'Units', 'Normalized');

    hold on;
    p1 = pcolor(mTime, height, visMat_Holger * 1e-3);
    p1.EdgeColor = 'none';

    xlabel('时间');
    ylabel('距离 (千米)');
    title(sprintf('%s 能见度 (前向迭代)', datestr(mTime(1), 'yyyy-mm-dd')));

    xlim([min(mTime), max(mTime)]);
    ylim([0, 10]);
    caxis([0, 50]);

    colormap('jet');

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'XTick', linspace(min(mTime), max(mTime), 5), 'layer', 'Top', 'TickDir', 'out', 'FontSize', 11);

    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    cb = colorbar('Position', [0.85, 0.3, 0.03, 0.5], 'Units', 'Normalized');
    titleHandle = get(cb, 'Title');
    set(titleHandle, 'String', '[千米]');
    set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

    export_fig(gcf, fullfile(saveFolder, sprintf('%s_vis_thi_Holger.png', datestr(mTime(1), 'yyyy-mm-dd'))), '-r300');

    %% Fernald
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w');

    subplot('Position', [0.1, 0.2, 0.72, 0.7], 'Units', 'Normalized');

    hold on;
    p1 = pcolor(mTime, height, visMat_Fernald * 1e-3);
    p1.EdgeColor = 'none';

    xlabel('时间');
    ylabel('距离 (千米)');
    title(sprintf('%s 能见度 (Fernald算法)', datestr(mTime(1), 'yyyy-mm-dd')));

    xlim([min(mTime), max(mTime)]);
    ylim([0, 10]);
    caxis([0, 50]);

    colormap('jet');

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'XTick', linspace(min(mTime), max(mTime), 5), 'layer', 'Top', 'TickDir', 'out', 'FontSize', 11);

    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    cb = colorbar('Position', [0.85, 0.3, 0.03, 0.5], 'Units', 'Normalized');
    titleHandle = get(cb, 'Title');
    set(titleHandle, 'String', '[千米]');
    set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

    export_fig(gcf, fullfile(saveFolder, sprintf('%s_vis_thi_Fernald.png', datestr(mTime(1), 'yyyy-mm-dd'))), '-r300');

    %% LK
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w');

    subplot('Position', [0.1, 0.2, 0.72, 0.7], 'Units', 'Normalized');

    hold on;
    p1 = pcolor(mTime, height, visMat_LK * 1e-3);
    p1.EdgeColor = 'none';

    xlabel('时间');
    ylabel('距离 (千米)');
    title(sprintf('%s 能见度 (LK算法)', datestr(mTime(1), 'yyyy-mm-dd')));

    xlim([min(mTime), max(mTime)]);
    ylim([0, 10]);
    caxis([0, 50]);

    colormap('jet');

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'XTick', linspace(min(mTime), max(mTime), 5), 'layer', 'Top', 'TickDir', 'out', 'FontSize', 11);

    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    cb = colorbar('Position', [0.85, 0.3, 0.03, 0.5], 'Units', 'Normalized');
    titleHandle = get(cb, 'Title');
    set(titleHandle, 'String', '[千米]');
    set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

    export_fig(gcf, fullfile(saveFolder, sprintf('%s_vis_thi_LK.png', datestr(mTime(1), 'yyyy-mm-dd'))), '-r300');

    %% Comparison

    %% I0297
    refIdx297 = find(height >= dist2_297, 1);
    visFernald297 = visMat_Fernald(refIdx297, :);
    visXian297 = visMat_Xian(refIdx297, :);
    visHolger297 = visMat_Holger(refIdx297, :);
    visLK297 = visMat_LK(refIdx297, :);
    vis297Interp = interp1(vis297.mTime, vis297.vis, mTime);

    diffFernald297 = visFernald297 - vis297Interp;
    stdFernald297 = nanstd(diffFernald297((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)));
    meanFernald297 = nanmean(diffFernald297((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)));
    relStdFernald297 = nanstd(diffFernald297((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)) ./ vis297Interp((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)));
    relMeanFernald297 = nanmean(diffFernald297((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)) ./ vis297Interp((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)));
    nFernald297Outliers = sum(abs(diffFernald297) > visDiffThresh);

    diffXian297 = visXian297 - vis297Interp;
    stdXian297 = nanstd(diffXian297((abs(diffXian297) <= visDiffThresh) & (visXian297 > 0)));
    meanXian297 = nanmean(diffXian297((abs(diffXian297) <= visDiffThresh) & (visXian297 > 0)));
    relStdXian297 = nanstd(diffXian297((abs(diffXian297) <= visDiffThresh) & (visXian297 > 0)) ./ vis297Interp((abs(diffXian297) <= visDiffThresh) & (visXian297 > 0)));
    relMeanXian297 = nanmean(diffXian297((abs(diffXian297) <= visDiffThresh) & (visXian297 > 0)) ./ vis297Interp((abs(diffXian297) <= visDiffThresh) & (visXian297 > 0)));
    nXian297Outliers = sum(abs(diffXian297) > visDiffThresh);

    diffHolger297 = visHolger297 - vis297Interp;
    stdHolger297 = nanstd(diffHolger297((abs(diffHolger297) <= visDiffThresh) & (visHolger297 > 0)));
    meanHolger297 = nanmean(diffHolger297((abs(diffHolger297) <= visDiffThresh) & (visHolger297 > 0)));
    relStdHolger297 = nanstd(diffHolger297((abs(diffHolger297) <= visDiffThresh) & (visHolger297 > 0)) ./ vis297Interp((abs(diffHolger297) <= visDiffThresh) & (visHolger297 > 0)));
    relMeanHolger297 = nanmean(diffHolger297((abs(diffHolger297) <= visDiffThresh) & (visHolger297 > 0)) ./ vis297Interp((abs(diffHolger297) <= visDiffThresh) & (visHolger297 > 0)));
    nHolger297Outliers = sum(abs(diffHolger297) > visDiffThresh);

    diffLK297 = visLK297 - vis297Interp;
    stdLK297 = nanstd(diffLK297((abs(diffLK297) <= visDiffThresh) & (visLK297 > 0)));
    meanLK297 = nanmean(diffLK297((abs(diffLK297) <= visDiffThresh) & (visLK297 > 0)));
    relStdLK297 = nanstd(diffLK297((abs(diffLK297) <= visDiffThresh) & (visLK297 > 0)) ./ vis297Interp((abs(diffLK297) <= visDiffThresh) & (visLK297 > 0)));
    relMeanLK297 = nanmean(diffLK297((abs(diffLK297) <= visDiffThresh) & (visLK297 > 0)) ./ vis297Interp((abs(diffLK297) <= visDiffThresh) & (visLK297 > 0)));
    nLK297Outliers = sum(abs(diffLK297) > visDiffThresh);

    % write statistics
    fprintf(sFid, 'Date: %s\n', datestr(mTime(1), 'yyyy-mm-dd'));
    fprintf(sFid, 'Total Profiles: %d\n\n', length(mTime));
    fprintf(sFid, 'I0297\n');
    fprintf(sFid, '方法: 平均偏差 平均相对偏差 标准差 平均相对偏差标准差 异常点个数(偏差超过%6.2fkm)\n', visDiffThresh * 1e-3);
    fprintf(sFid, 'Fernald: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanFernald297 * 1e-3, relMeanFernald297 * 100, stdFernald297 * 1e-3, relStdFernald297 * 100, nFernald297Outliers);
    fprintf(sFid, 'Xian: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanXian297 * 1e-3, relMeanXian297 * 100, stdXian297 * 1e-3, relStdXian297 * 100, nXian297Outliers);
    fprintf(sFid, 'Holger: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanHolger297 * 1e-3, relMeanHolger297 * 100, stdHolger297 * 1e-3, relStdHolger297 * 100, nHolger297Outliers);
    fprintf(sFid, 'LK: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanLK297 * 1e-3, relMeanLK297 * 100, stdLK297 * 1e-3, relStdLK297 * 100, nLK297Outliers);

    figure('Position', [0, 30, 550, 500], 'Units', 'Pixels', 'Color', 'w');

    subplot('Position', [0.15, 0.6, 0.8, 0.35], 'Units', 'Normalized');

    hold on;
    p1 = plot(mTime, visMat_Fernald(refIdx297, :) * 1e-3, '-k', 'DisplayName', 'Fernald算法');
    p2 = plot(mTime, visMat_Xian(refIdx297, :) * 1e-3, '-b', 'DisplayName', '大舜算法');
    p3 = plot(mTime, visMat_Holger(refIdx297, :) * 1e-3, '-g', 'DisplayName', '前向迭代');
    p4 = plot(mTime, visMat_LK(refIdx297, :) * 1e-3, '-r', 'DisplayName', '蓝科光电');
    p5 = scatter(vis297.mTime, vis297.vis * 1e-3, 10, 'Marker', 's', 'MarkerEdgeColor', 'cyan', 'MarkerFaceColor', 'cyan', 'DisplayName', 'I0297');
    hold off;

    xlim([min(mTime), max(mTime)]);
    ylim([0, 50]);

    % xlabel('时间');
    ylabel('能见度 (千米)');
    title(sprintf('%s 能见度对比 (I0297)', datestr(mTime(1), 'yyyy-mm-dd')));

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    legend([p1, p2, p3, p4, p5], 'Location', 'NorthEast');

    subplot('Position', [0.15, 0.1, 0.8, 0.25], 'Units', 'Normalized');

    hold on;

    plot(mTime, diffFernald297 * 1e-3, '-k');
    plot(mTime, diffXian297 * 1e-3, '-b');
    plot(mTime, diffHolger297 * 1e-3, '-g');
    plot(mTime, diffLK297 * 1e-3, '-r');

    plot(mTime, zeros(size(mTime)), 'LineStyle', '-.', 'color', 'cyan');

    hold off;

    xlim([min(mTime), max(mTime)]);
    ylim([-10, 10]);

    xlabel('时间');
    ylabel('能见度偏差 (千米)');

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    text(0, 1.24, sprintf('(Fernald 算法) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(大舜算法) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(前向迭代) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(蓝科光电) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n', ...
        meanFernald297 * 1e-3, relMeanFernald297 * 100, stdFernald297 * 1e-3, relStdFernald297 * 100, ...
        meanXian297 * 1e-3, relMeanXian297 * 100, stdXian297 * 1e-3, relStdXian297 * 100, ...
        meanHolger297 * 1e-3, relMeanHolger297 * 100, stdHolger297 * 1e-3, relStdHolger297 * 100, ...
        meanLK297 * 1e-3, relMeanLK297 * 100, stdLK297 * 1e-3, relStdLK297 * 100), 'Units', 'Normalized', 'FontSize', 10);

    export_fig(gcf, fullfile(saveFolder, sprintf('%s_vis_comparison_297.png', datestr(mTime(1), 'yyyy-mm-dd'))), '-r300');

    % laser temperature
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w');

    hold on;
    plot(mTime, laserTemp, '-k');
    hold off;

    xlim([min(mTime), max(mTime)]);
    ylim([0, 50]);

    xlabel('时间');
    ylabel('温度（度）');
    title(sprintf('%s 激光头温度', datestr(mTime(1), 'yyyy-mm-dd')));

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    export_fig(gcf, fullfile(saveFolder, sprintf('%s_laser_temperature.png', datestr(mTime(1), 'yyyy-mm-dd'))), '-r300');

    %% I0998
    refIdx998 = find(height >= dist2_998, 1);
    visFernald998 = visMat_Fernald(refIdx998, :);
    visXian998 = visMat_Xian(refIdx998, :);
    visHolger998 = visMat_Holger(refIdx998, :);
    visLK998 = visMat_LK(refIdx998, :);
    vis998Interp = interp1(vis998.mTime, vis998.vis, mTime);

    diffFernald998 = visFernald998 - vis998Interp;
    stdFernald998 = nanstd(diffFernald998((abs(diffFernald998) <= visDiffThresh) & (visFernald998 > 0)));
    meanFernald998 = nanmean(diffFernald998((abs(diffFernald998) <= visDiffThresh) & (visFernald998 > 0)));
    relStdFernald998 = nanstd(diffFernald998((abs(diffFernald998) <= visDiffThresh) & (visFernald998 > 0)) ./ vis998Interp((abs(diffFernald998) <= visDiffThresh) & (visFernald998 > 0)));
    relMeanFernald998 = nanmean(diffFernald998((abs(diffFernald998) <= visDiffThresh) & (visFernald998 > 0)) ./ vis998Interp((abs(diffFernald998) <= visDiffThresh) & (visFernald998 > 0)));
    nFernald998Outliers = sum(abs(diffFernald998) > visDiffThresh);

    diffXian998 = visXian998 - vis998Interp;
    stdXian998 = nanstd(diffXian998((abs(diffXian998) <= visDiffThresh) & (visXian998 > 0)));
    meanXian998 = nanmean(diffXian998((abs(diffXian998) <= visDiffThresh) & (visXian998 > 0)));
    relStdXian998 = nanstd(diffXian998((abs(diffXian998) <= visDiffThresh) & (visXian998 > 0)) ./ vis998Interp((abs(diffXian998) <= visDiffThresh) & (visXian998 > 0)));
    relMeanXian998 = nanmean(diffXian998((abs(diffXian998) <= visDiffThresh) & (visXian998 > 0)) ./ vis998Interp((abs(diffXian998) <= visDiffThresh) & (visXian998 > 0)));
    nXian998Outliers = sum(abs(diffXian998) > visDiffThresh);

    diffHolger998 = visHolger998 - vis998Interp;
    stdHolger998 = nanstd(diffHolger998((abs(diffHolger998) <= visDiffThresh) & (visHolger998 > 0)));
    meanHolger998 = nanmean(diffHolger998((abs(diffHolger998) <= visDiffThresh) & (visHolger998 > 0)));
    relStdHolger998 = nanstd(diffHolger998((abs(diffHolger998) <= visDiffThresh) & (visHolger998 > 0)) ./ vis998Interp((abs(diffHolger998) <= visDiffThresh) & (visHolger998 > 0)));
    relMeanHolger998 = nanmean(diffHolger998((abs(diffHolger998) <= visDiffThresh) & (visHolger998 > 0)) ./ vis998Interp((abs(diffHolger998) <= visDiffThresh) & (visHolger998 > 0)));
    nHolger998Outliers = sum(abs(diffHolger998) > visDiffThresh);

    diffLK998 = visLK998 - vis998Interp;
    stdLK998 = nanstd(diffLK998((abs(diffLK998) <= visDiffThresh) & (visLK998 > 0)));
    meanLK998 = nanmean(diffLK998((abs(diffLK998) <= visDiffThresh) & (visLK998 > 0)));
    relStdLK998 = nanstd(diffLK998((abs(diffLK998) <= visDiffThresh) & (visLK998 > 0)) ./ vis998Interp((abs(diffLK998) <= visDiffThresh) & (visLK998 > 0)));
    relMeanLK998 = nanmean(diffLK998((abs(diffLK998) <= visDiffThresh) & (visLK998 > 0)) ./ vis998Interp((abs(diffLK998) <= visDiffThresh) & (visLK998 > 0)));
    nLK998Outliers = sum(abs(diffLK998) > visDiffThresh);

    % write statistics
    fprintf(sFid, 'I0998\n');
    fprintf(sFid, '方法: 平均偏差 平均相对偏差 标准差 平均相对偏差标准差 异常点个数(偏差超过%6.2fkm)\n', visDiffThresh * 1e-3);
    fprintf(sFid, 'Fernald: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanFernald998 * 1e-3, relMeanFernald998 * 100 , stdFernald998 * 1e-3, relStdFernald998 * 100, nFernald998Outliers);
    fprintf(sFid, 'Xian: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanXian998 * 1e-3, relMeanXian998 * 100, stdXian998 * 1e-3, relStdXian998 * 100, nXian998Outliers);
    fprintf(sFid, 'Holger: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanHolger998 * 1e-3, relMeanHolger998 * 100, stdHolger998 * 1e-3, relStdHolger998 * 100, nHolger998Outliers);
    fprintf(sFid, 'LK: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n\n', meanLK998 * 1e-3, relMeanLK998 * 100, stdLK998 * 1e-3, relStdLK998 * 100, nLK998Outliers);

    figure('Position', [0, 30, 550, 500], 'Units', 'Pixels', 'Color', 'w');

    subplot('Position', [0.15, 0.6, 0.8, 0.35], 'Units', 'Normalized');

    hold on;
    p1 = plot(mTime, visMat_Fernald(refIdx998, :) * 1e-3, '-k', 'DisplayName', 'Fernald算法');
    p2 = plot(mTime, visMat_Xian(refIdx998, :) * 1e-3, '-b', 'DisplayName', '大舜算法');
    p3 = plot(mTime, visMat_Holger(refIdx998, :) * 1e-3, '-g', 'DisplayName', '前向迭代');
    p4 = plot(mTime, visMat_LK(refIdx998, :) * 1e-3, '-r', 'DisplayName', '蓝科光电');
    p5 = scatter(vis998.mTime, vis998.vis * 1e-3, 10, 'Marker', 's', 'MarkerEdgeColor', 'cyan', 'MarkerFaceColor', 'cyan', 'DisplayName', 'I0998');
    hold off;

    xlim([min(mTime), max(mTime)]);
    ylim([0, 50]);

    % xlabel('时间');
    ylabel('能见度 (千米)');
    title(sprintf('%s 能见度对比 (I0998)', datestr(mTime(1), 'yyyy-mm-dd')));

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    legend([p1, p2, p3, p4, p5], 'Location', 'NorthEast');

    subplot('Position', [0.15, 0.1, 0.8, 0.25], 'Units', 'Normalized');

    hold on;

    plot(mTime, diffFernald998 * 1e-3, '-k');
    plot(mTime, diffXian998 * 1e-3, '-b');
    plot(mTime, diffHolger998 * 1e-3, '-g');
    plot(mTime, diffLK998 * 1e-3, '-r');

    plot(mTime, zeros(size(mTime)), 'LineStyle', '-.', 'color', 'cyan');

    hold off;

    xlim([min(mTime), max(mTime)]);
    ylim([-10, 10]);

    xlabel('时间');
    ylabel('能见度偏差 (千米)');

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    text(0, 1.24, sprintf('(Fernald 算法) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(大舜算法) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(前向迭代) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(蓝科光电) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n', ...
        meanFernald998 * 1e-3, relMeanFernald998 * 100, stdFernald998 * 1e-3, relStdFernald998 * 100, ...
        meanXian998 * 1e-3, relMeanXian998 * 100, stdXian998 * 1e-3, relStdXian998 * 100, ...
        meanHolger998 * 1e-3, relMeanHolger998 * 100, stdHolger998 * 1e-3, relStdHolger998 * 100, ...
        meanLK998 * 1e-3, relMeanLK998 * 100, stdLK998 * 1e-3, relStdLK998 * 100), 'Units', 'Normalized', 'FontSize', 10);

    export_fig(gcf, fullfile(saveFolder, sprintf('%s_vis_comparison_998.png', datestr(mTime(1), 'yyyy-mm-dd'))), '-r300');

end

fclose(sFid);