%% Parameter Definition
visSensorFile = 'vis-sensor-data.mat';   % 前向散射能见度仪数据文件
l0Folder = 'G:\backup\vis-lidar\20260428-tianjin-evaluation\前散仪对比数据\VIS1\L0';   % 雷达原始数据文件目录
l1Folder = 'G:\backup\vis-lidar\20260428-tianjin-evaluation\前散仪对比数据\VIS1\L1';   % 雷达产品文件目录
saveFolder = 'H:\research\vislidar-intercomparison\tianjing';   % 输出结果目录
overlapFile = 'VIS1_overlap_factor.txt';   % 重叠因子文件，重叠因子通过水平方法计算
flagOLCor = true;   % 是否进行重叠因子修正
zenithAgl = 6;   % 雷达仰角（度）
lr = 50;   % 雷达比
refH = [8, 9];   % Fernald方法参考高度（千米）
flagReadData = true;
AEConvFactor = (1030/550) ^1;   % 波长指数（转换因子）
flagDisplay = 'on';

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

    thisDate = datenum(basename(dateFolders{iFolder}), 'yyyymmdd');

    % load level 0 data
    l0DataFile = fullfile(l0Folder, sprintf('%s_vis_lidar_l0.mat', datestr(thisDate, 'yyyy-mm-dd')));
    load(l0DataFile);

    % % load level 1 data
    % l1DataFile = fullfile(l1Folder, sprintf('%s_vis_lidar_l1.mat', datestr(thisDate, 'yyyy-mm-dd')));
    % load(l1DataFile);

    % preprocessing
    thisSig = lidarSig(:, 1);
    thisBG = nanmean(thisSig((end - 200):end));
    thisSigCor = thisSig - thisBG;
    if flagOLCor
        ol = interp1(olHeight, olVal, range);
        thisSigCor = thisSigCor ./ ol;
    end
    thisRCS = thisSigCor .* range.^2;

    % Rayleigh scattering
    [mBsc, mExt] = MolModel(range * sin(zenithAgl / 180 * pi), 1030, 'meteor', 'standard_atmosphere');
    [~, mExt550] = MolModel(range * sin(zenithAgl / 180 * pi), 550, 'meteor', 'standard_atmosphere');
    mAttn = mBsc .* exp(-2 * nancumsum(mExt .* [range(1); diff(range)] * 1e3));
    ratioL2M = nansum(thisRCS(1500:1550)) / nansum(mAttn(1500:1550));

    % lidar calibration
    isInHCaliRange = (range * 1e-3 >= refH(1)) & (range * 1e-3 <= refH(2));
    try
        [~, extRef] = chi2fit(range(isInHCaliRange) * 1e3, log(thisRCS(isInHCaliRange)), 0.1 * ones(size(thisRCS(isInHCaliRange))));
    catch
        extRef = 1e-5;
    end
    extRef = extRef / -2;
    aBsc = transpose(fernald(range * 1e3, thisSigCor, thisBG, lr, refH * 1e3, extRef / lr, mBsc, 4));

    aBscCor = aBsc;
    aBscCor(range <= 1.5) = aBscCor(find(range >= 1.5, 1));
    lc = thisRCS ./ ((aBscCor + mBsc) .* exp(-2 * nancumsum(aBscCor * lr + mExt).* [range(1); diff(range)] * 1e3));
    lcMean = nanmean(lc(500:600));

    %% signal preprocessing
    thisBG = nanmean(lidarSig((end - 200):end, :), 1);
    lidarSigCor = lidarSig - repmat(thisBG, size(thisBG, 1), 1);
    if flagOLCor
        lidarSigCor = lidarSigCor ./ repmat(ol, 1, size(lidarSigCor, 2));
    end
    lidarRCS = lidarSigCor .* repmat(range, 1, size(lidarSigCor, 2)).^2;

    %% Display

    % rcs
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w', 'visible', flagDisplay);

    hold on;
    p1 = pcolor(mTime, range * 1e-3, lidarRCS);
    p1.EdgeColor = 'none';
    hold off;

    xlabel('时间');
    ylabel('距离 (千米)');
    title('距离修正信号');

    xlim([min(mTime), max(mTime)]);
    ylim([0, 10]);
    caxis([0, 5e10]);

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'layer', 'Top', 'TickDir', 'out', 'FontSize', 11);

    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    colorbar();
    export_fig(gcf, fullfile(saveFolder, sprintf('%s_RCS.png', datestr(mTime(1), 'yyyy-mm-dd'))), '-r300');

    %% extinction
    extMat_Xian = NaN(size(lidarSigCor));
    extMat_Holger = NaN(size(lidarSigCor));
    extMat_Fernald = NaN(size(lidarSigCor));
    [MTIME, HEIGHT] = meshgrid(mTime, range);
    for iPrf = 1:size(lidarSigCor, 2)
        thisRCS = lidarRCS(:, iPrf);
        thisRCS(thisRCS <= 0) = NaN;
        isInHCaliRange = (range >= refH(1)) & (range <= refH(2));
        [~, extRef] = chi2fit(range(isInHCaliRange) * 1e3, log(thisRCS(isInHCaliRange)), 0.1 * ones(size(thisRCS(isInHCaliRange))));
        extRef = extRef / -2;

        % extinction retrieval
        ext1 = extRet_Fernald(range, lidarSigCor(:, iPrf), thisBG(iPrf), mBsc, mExt, 'minSNR', 3, 'hFullOL', 330, 'calibration_constant', lcMean * 1e6);
        [~, ext2] = extRet_Holger(range', lidarSigCor(:, iPrf)', ...
                                'calibration_constant', lcMean * 1e6, ...
                                'fullOverlapR', 1200, ...
                                'elevation_angle', 0);
        ext3 = extRet_Xian(range, lidarSigCor(:, iPrf), thisBG(iPrf), 'minSNR', 2, 'rangeFullOverlap', 330);

        extMat_Xian(:, iPrf) = ext3;
        extMat_Holger(:, iPrf) = ext2;
        extMat_Fernald(:, iPrf) = ext1;
    end

    %% extinction to visibility
    visMat_Xian = ext2vis(extMat_Xian);
    visMat_Holger = ext2vis(extMat_Holger * AEConvFactor + mExt550);
    visMat_Fernald = ext2vis(extMat_Fernald * AEConvFactor + mExt550);

    %% Save to .mat file
    save(fullfile(l1Folder, sprintf('%s_vis_lidar_l1_exp.mat', datestr(mTime(1), 'yyyy-mm-dd'))), ...
        'mTime', 'range', 'visMat_Xian', 'visMat_Holger', 'visMat_Fernald', ...
        'extMat_Xian', 'extMat_Holger', 'extMat_Fernald');

    %% Display

    % Xian
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w', 'visible', flagDisplay);

    subplot('Position', [0.1, 0.2, 0.72, 0.7], 'Units', 'Normalized');

    hold on;
    p1 = pcolor(mTime, range * 1e-3, visMat_Xian * 1e-3);
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
    p1 = pcolor(mTime, range * 1e-3, visMat_Holger * 1e-3);
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
    p1 = pcolor(mTime, range * 1e-3, visMat_Fernald * 1e-3);
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

end