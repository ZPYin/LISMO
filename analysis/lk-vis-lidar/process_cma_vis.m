%% Parameter Definition
visSensorFile = 'vis-sensor-data.mat';   % 前向散射能见度仪数据文件
l0Folder = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\L0';   % 雷达原始数据文件目录
l1Folder = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\L1';   % 雷达产品文件目录
saveFolder = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\quicklooks';   % 输出结果目录
overlapFile = 'lk_20241101_overlap_factor.txt';   % 重叠因子文件，重叠因子通过水平方法计算
flagOLCor = true;   % 是否进行重叠因子修正
lr = 50;   % 雷达比
dist2_998 = 4.5;   % I0998站点与激光雷达距离（千米）
dist2_297 = 3.2;   % I0297站点与激光雷达距离（千米）
refH = [8, 9];   % Fernald方法参考高度（千米）
nPretrigger = 12;
flagReadData = true;
hRes = 7.5;   % 高度分辨率（米）
% AEConvFactor = (1064/550) ^1.0915;   % 波长指数（转换因子）
AEConvFactor = (1064/550) ^1;   % 波长指数（转换因子）
visDiffThresh = 2e4;   % [m]
flagDisplay = 'off';

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


%% find dates
dateFolders = listdir(l0Folder, '\w*', 1);

for iFolder = 1:length(dateFolders)

    fprintf('Finished %6.2f%%: processing %s\n', (iFolder - 1) / length(dateFolders) * 100, dateFolders{iFolder});

    close all;

    thisDate = datenum(basename(dateFolders{iFolder}), 'yyyy-mm-dd');

    % load level 0 data
    l0DataFile = fullfile(l0Folder, sprintf('%s_vis_lidar_l0.mat', datestr(thisDate, 'yyyy-mm-dd')));
    load(l0DataFile);

    % load level 1 data
    l1DataFile = fullfile(l1Folder, sprintf('%s_vis_lidar_l1.mat', datestr(thisDate, 'yyyy-mm-dd')));
    load(l1DataFile);

    % preprocessing
    thisSig = lidarSig(:, 1);
    thisBG = nanmean(thisSig((end - 200):end));
    thisSigCor = thisSig - thisBG;
    height = transpose(((1:length(thisSigCor)) - 0.5) * hRes * 1e-3);
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

    %% signal preprocessing
    thisBG = nanmean(lidarSig((end - 200):end, :), 1);
    lidarSigCor = lidarSig - repmat(thisBG, size(thisBG, 1), 1);
    if flagOLCor
        lidarSigCor = lidarSigCor ./ repmat(ol, 1, size(lidarSigCor, 2));
    end
    lidarRCS = lidarSigCor .* repmat(height, 1, size(lidarSigCor, 2)).^2;

    %% Display

    % rcs
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w', 'visible', flagDisplay);

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

    %% extinction
    extMat_Xian = NaN(size(lidarSigCor));
    extMat_Holger = NaN(size(lidarSigCor));
    extMat_Fernald = NaN(size(lidarSigCor));
    [MTIMELK, HEIGHTLK] = meshgrid(mTime_LK, height_LK * 1e-3);
    [MTIME, HEIGHT] = meshgrid(mTime, height);
    extMat_LK = interp2(mTime_LK, height_LK, ext_LK, MTIME, HEIGHT);
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

    %% extinction to visibility
    visMat_Xian = ext2vis(extMat_Xian);
    visMat_Holger = ext2vis(extMat_Holger * AEConvFactor + mExt550);
    visMat_Fernald = ext2vis(extMat_Fernald * AEConvFactor + mExt550);
    visMat_LK = interp2(MTIMELK, HEIGHTLK, vis_LK * 1e3, MTIME, HEIGHT);

    %% Save to .mat file
    save(fullfile(l1Folder, sprintf('%s_vis_lidar_l1_exp.mat', datestr(mTime(1), 'yyyy-mm-dd'))), ...
        'mTime', 'height', 'visMat_Xian', 'visMat_Holger', 'visMat_Fernald', 'visMat_LK', ...
        'extMat_Xian', 'extMat_Holger', 'extMat_Fernald', 'extMat_LK');

    %% Display

    % Xian
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w', 'visible', flagDisplay);

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
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w', 'visible', flagDisplay);

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
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w', 'visible', flagDisplay);

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
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w', 'visible', flagDisplay);

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

end