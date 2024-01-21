clc;
close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\RawData\QingDao1';
scanTime2Plot = datenum(2023, 12, 27, 22, 23, 0);   % temporal range for data input
olFile = 'overlap_20231227.mat';
distOffset = -48.75;
iCh = 1;
flagReadScanData = true;
flagReadAllDayData = true;
flagAve4Vis = false;
iPrfInScan4Vis = 1;
flagOverlapCor = true;
olHeight = 1000;
refDist = [3000, 4500];
calDist = [2000, 3000];
visRetMethod = 'xian';
overlapFile = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\wxzk\near-sea-shore-measurement\2023-12-27\overlap_20231227.mat';
savePath = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\wxzk\near-sea-shore-measurement\2023-12-27';
visFile = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\wxzk\near-sea-shore-measurement\vis-fixed-platform.mat';

%% Read Data
dayPath = fullfile(dataFolder, datestr(scanTime2Plot(1), 'yyyy'), datestr(scanTime2Plot(1), 'mm'), datestr(scanTime2Plot(1), 'dd'));
scanPaths = listdir(dayPath, '.*', 1);

% Read Scan Data
if flagReadScanData

    scanTimes = [];
    for iScanPath = 1:length(scanPaths)
        basePathStr = basename(scanPaths{iScanPath});

        % H_SCAN15_135_2_20231226001201
        basePathTime = datenum(basePathStr((end - 13):end), 'yyyymmddHHMMSS');
        scanTimes = cat(2, scanTimes, basePathTime);
    end
    [~, iPath] = min(abs(scanTimes - scanTime2Plot));

    scanData = readVIS(scanPaths{iPath}, 'isDir', true, 'debug', false);
end

% Read Full Day Data
if flagReadAllDayData
    fullData = struct();
    fullData.mTime = [];
    fullData.rawSignal = [];
    fullData.nShots = [];

    for iScan = 1:length(scanPaths)
        fprintf('Finished %6.2f%%: reading %s\n', (iScan - 1) / length(scanPaths) * 100, scanPaths{iScan});

        scanFiles = listfile(scanPaths{iScan}, '\w*', 1);
        if ~ flagAve4Vis
            scanFiles = scanFiles{iPrfInScan4Vis};
        end

        lData = readVIS(scanFiles, 'isDir', false);

        fullData.mTime = cat(2, fullData.mTime, mean(lData.startTime));
        fullData.nShots = cat(2, fullData.nShots, sum(lData.nShots));
        fullData.rawSignal = cat(1, fullData.rawSignal, sum(lData.rawSignal, 1));
    end
end

% Read Vis Data
visData = load(visFile);

%% Preprocess
range = ((1:lData.nBins(1)) - 0.5) * lData.hRes(1) + distOffset;

if exist('scanData', 'var')
    scanData.bg = nanmean(squeeze(scanData.rawSignal(:, iCh, (end - 20):end)), 2);
    scanData.signal = squeeze(scanData.rawSignal(:, iCh, :)) - repmat(scanData.bg, 1, scanData.nBins(1));
    scanData.rcs = scanData.signal .* repmat(range, size(scanData.rawSignal, 1), 1).^2;
    scanData.snr = scanData.signal ./ sqrt(squeeze(scanData.rawSignal(:, iCh, :)));
    scanData.lowSNRMask = false(size(scanData.signal));
    for iPrf = 1:size(scanData.signal, 1)
        snrPrf = scanData.snr(iPrf, :);
        snrPrf(range < olHeight) = NaN;
        snrLow = find(snrPrf < 1, 1);
        scanData.lowSNRMask(iPrf, snrLow:end) = true;
    end
end

if exist('fullData', 'var')
    fullData.bg = nanmean(squeeze(fullData.rawSignal(:, iCh, (end - 20):end)), 2);
    fullData.signal = squeeze(fullData.rawSignal(:, iCh, :)) - repmat(fullData.bg, 1, scanData.nBins(1));
    fullData.rcs = fullData.signal .* repmat(range, size(fullData.rawSignal, 1), 1).^2;
    fullData.snr = fullData.signal ./ sqrt(squeeze(fullData.rawSignal(:, iCh, :)));
    fullData.lowSNRMask = false(size(fullData.signal));
    for iPrf = 1:size(fullData.signal, 1)
        snrPrf = fullData.snr(iPrf, :);
        snrPrf(range < olHeight) = NaN;
        snrLow = find(snrPrf < 1, 1);
        fullData.lowSNRMask(iPrf, snrLow:end) = true;
    end
end

%% Overlap Correction
if flagOverlapCor
    ol = load(overlapFile);

    if exist('fullData', 'var')
        fullData.signal = fullData.signal ./ repmat(transpose(ol.ov), size(fullData.signal, 1), 1);
        fullData.rcs = fullData.rcs ./ repmat(transpose(ol.ov), size(fullData.rcs, 1), 1);
    end
    
    if exist('scanData', 'var')
        scanData.signal = scanData.signal ./ repmat(transpose(ol.ov), size(scanData.signal, 1), 1);
        scanData.rcs = scanData.rcs ./ repmat(transpose(ol.ov), size(scanData.rcs, 1), 1);
    end
end

%% lidar calibration
meanRCS = nanmean(scanData.rcs, 1);
height = range .* cos(lData.zenithAng(1) / 180 * pi);
[mBsc, mExt] = MolModel(height, 1064, 'meteor', 'standard_atmosphere');   % Rayleigh Scattering
mAttn = mBsc .* exp(-2 * nancumsum(mExt .* [range(1), diff(range)]));
isInRefDist = (range >= refDist(1)) & (range <= refDist(2));
ratioL2M = nanmean(mAttn(isInRefDist)) / nanmean(meanRCS(isInRefDist));

% calculate reference value
extSlopeMethod = movingslope(log(smooth(meanRCS, 8)), 8) / (-2) / (range(2) - range(1)) - transpose(mExt);
extRef = nanmean(extSlopeMethod(isInRefDist));
extRefStd = nanstd(extSlopeMethod(isInRefDist));
aBsc = fernald(range, nanmean(scanData.signal, 1), nansum(scanData.bg), 50, refDist, extRef / 50, mBsc, 4);
aBscCor = aBsc;
aBscCor(range <= olHeight) = aBscCor(find(range >= olHeight, 1));
lc = meanRCS ./ ((aBscCor + mBsc) .* exp(-2 * nancumsum((aBscCor .* 50 + mExt) .* [range(1), diff(range)])));
lcMean = nanmean(lc(isInRefDist));
lcStd = nanstd(lc(isInRefDist));
meanL2M = nanmean(meanRCS(isInRefDist)) / nanmean(mAttn(isInRefDist));

%% extinction&visibility retrieval
if exist('scanData', 'var')
    scanData.ext = NaN(size(scanData.signal));
    scanData.vis = NaN(size(scanData.signal));

    for iPrf = 1:size(scanData.signal, 1)

        if strcmpi(visRetMethod, 'xian')
            scanData.ext(iPrf, :) = extRet_Xian(range, scanData.signal(iPrf, :), scanData.bg(iPrf), 'minSNR', 0.5, 'rangeFullOverlap', 200);
        elseif strcmpi(visRetMethod, 'quasi')
            [~, scanData.ext(iPrf, :)] = extRet_Holger(range, scanData.signal(iPrf, :), ...
                'calibration_constant', lcMean, ...
                'fullOverlapR', 100, ...
                'elevation_angle', scanData.zenithAng(iPrf));
            scanData.ext(iPrf, :) = scanData.ext(iPrf, :) + mExt;
        else
            warning('Undefined extinction retrieval method.');
        end

    end

    scanData.ext(scanData.lowSNRMask) = NaN;
    scanData.vis = ext2vis(scanData.ext);
    scanData.vis(isnan(scanData.vis)) = 1e5;
end

if exist('fullData', 'var')
    fullData.ext = NaN(size(fullData.signal));
    fullData.vis = NaN(size(fullData.signal));

    for iPrf = 1:size(fullData.signal, 1)

        if strcmpi(visRetMethod, 'xian')
            fullData.ext(iPrf, :) = extRet_Xian(range, fullData.signal(iPrf, :), fullData.bg(iPrf), 'minSNR', 0.5, 'rangeFullOverlap', 200);
        elseif strcmpi(visRetMethod, 'quasi')
            [~, fullData.ext(iPrf, :)] = extRet_Holger(range, fullData.signal(iPrf, :), ...
                'calibration_constant', lcMean, ...
                'fullOverlapR', 100, ...
                'elevation_angle', scanData.zenithAng(1));
            fullData.ext(iPrf, :) = fullData.ext(iPrf, :) + mExt;
        else
            warning('Undefined extinction retrieval method.');
        end

    end

    fullData.ext(fullData.lowSNRMask) = NaN;
    fullData.vis = ext2vis(fullData.ext);
    fullData.vis(isnan(fullData.vis)) = 1e5;
end

%% Display

% scan rcs
if exist('scanData', 'var')
    figure('Position', [0, 0, 600, 400], 'color', 'w', 'visible', 'on');

    subplot('Position', [0.1, 0, 0.86, 0.86], 'Units', 'normalized');

    scanData.rcs(scanData.lowSNRMask) = NaN;
    [~, p1] = polarPcolor(range / 1e3, scanData.azimuthAng, scanData.rcs, 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 10], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 8e9], 'tickSize', 12, 'tickColor', 'm');
    ylabel(p1, 'range cor. sig.');
    set(p1, 'location', 'westoutside', 'FontSize', 12);
    colormap(gca, myColormap('jetImage'));
    text(0.3, 1.2, sprintf('%s', datestr(mean(scanData.startTime), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'Bold');
    text(0.6, 0, 'distance (km)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'Bold');
end

% calibration
if exist('scanData', 'var')
    figure('Position', [0, 20, 500, 500], 'Units', 'Pixels', 'Color', 'w');

    subfig = subfigPos([0.12, 0.1, 0.85, 0.8], 2, 1, 0, 0.1);

    subplot('Position', subfig(1, :), 'Units', 'normalized');

    hold on;

    p1 = plot(range * 1e-3, mAttn * 1e3, 'Color', [242, 89, 75] / 255, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', '分子信号');
    p2 = plot(range * 1e-3, meanRCS * ratioL2M * 1e3, 'Color', [65, 54, 89] / 255, 'LineWidth', 2, 'DisplayName', '雷达');

    plot([refDist(1), refDist(1)] * 1e-3, [1e-20, 1e20], '--k');
    plot([refDist(2), refDist(2)] * 1e-3, [1e-20, 1e20], '--k');

    xlabel('水平距离 (千米)');
    ylabel('衰减后向散射系数 (km-1sr-1)');

    xlim([0, 5]);
    ylim([1e-6, 1e-3]);

    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log');

    legend([p1, p2], 'Location', 'SouthEast');

    subplot('Position', subfig(2, :), 'Units', 'normalized');

    hold on;

    p1 = plot(range * 1e-3, mBsc * 1e3, 'Color', [242, 89, 75] / 255, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', '分子信号');
    p2 = plot(range * 1e-3, aBsc * 1e3 * 50, 'Color', [65, 54, 89] / 255, 'LineWidth', 2, 'DisplayName', '雷达');

    plot([refDist(1), refDist(1)] * 1e-3, [1e-20, 1e20], '--k');
    plot([refDist(2), refDist(2)] * 1e-3, [1e-20, 1e20], '--k');

    xlabel('水平距离 (千米)');
    ylabel('消光系数 (km-1)');

    xlim([0, 5]);
    ylim([-0.001, 1]);

    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

    figure('Position', [0, 30, 500, 300], 'Units', 'Pixels', 'Color', 'w');

    hold on;
    plot(range * 1e-3, lc, '-k', 'LineWidth', 2);
    plot([calDist(1), calDist(1)] * 1e-3, [1e-20, 1e20], '--k');
    plot([calDist(2), calDist(2)] * 1e-3, [1e-20, 1e20], '--k');

    text(0.3, 0.8, sprintf('雷达标定常数平均: %5.3e\\pm%5.3e\n参考分子标定常数:%5.3e\n', lcMean, lcStd, meanL2M), 'Units', 'normalized');

    xlabel('水平距离 (千米)');
    ylabel('雷达标定常数 (每单廓线)');

    xlim([0, 5]);
    ylim([0., lcMean * 3]);

    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');
end

% scan ext
if exist('scanData', 'var')
    figure('Position', [0, 0, 600, 400], 'color', 'w', 'visible', 'on');

    subplot('Position', [0.1, 0, 0.86, 0.86], 'Units', 'normalized');

    scanData.ext(scanData.lowSNRMask) = NaN;
    [~, p1] = polarPcolor(range / 1e3, scanData.azimuthAng, scanData.ext * 1e3, 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 10], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 10], 'tickSize', 12, 'tickColor', 'm');
    ylabel(p1, 'extinction (km-1)');
    set(p1, 'location', 'westoutside', 'FontSize', 12);
    colormap(gca, myColormap('jetImage'));
    text(0.3, 1.2, sprintf('%s', datestr(mean(scanData.startTime), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'Bold');
    text(0.6, 0, 'distance (km)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'Bold');
end

% scan vis
if exist('scanData', 'var')
    figure('Position', [0, 0, 600, 400], 'color', 'w', 'visible', 'on');

    subplot('Position', [0.1, 0, 0.86, 0.86], 'Units', 'normalized');

    scanData.vis(scanData.lowSNRMask) = NaN;
    [~, p1] = polarPcolor(range / 1e3, scanData.azimuthAng, scanData.vis * 1e-3, 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 10], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 30], 'tickSize', 12, 'tickColor', 'm');
    ylabel(p1, 'visibility (km)');
    set(p1, 'location', 'westoutside', 'FontSize', 12);
    load('vis_colormap.mat');
    colormap(gca, double(visColorbar) / 255);
    text(0.3, 1.2, sprintf('%s', datestr(mean(scanData.startTime), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'Bold');
    text(0.6, 0, 'distance (km)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'Bold');
end

% scan vis comparison
if exist('scanData', 'var')
    figure('Position', [0, 10, 550, 300], 'Units', 'Pixels', 'Color', 'w');

    hold on;
    s1 = plot(visData.mTime, visData.vis * 1e-3, 'color', 'k', 'Marker', '.', 'MarkerFaceColor', 'k', 'markeredgecolor', 'k', 'DisplayName', '能见度仪');
    s2 = plot(scanData.startTime, scanData.vis(:, 80) * 1e-3, 'color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b', 'markeredgecolor', 'b', 'DisplayName', sprintf('激光测雾雷达(%3d米)', floor(range(80))));
    hold off;

    xlim([scanData.startTime(1), scanData.startTime(end)]);
    ylim([0, 50]);

    xlabel('时间');
    ylabel('能见度 (千米)');
    title(sprintf('能见度序列对比(%s)', datestr(scanData.startTime(1), 'yyyy-mm-dd')));

    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);
    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');
    legend([s1, s2], 'location', 'NorthEast');
end

% full data rcs
if exist('fullData', 'var')

    figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

    subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
    hold on;
    rcsTmp = fullData.rcs;
    rcsTmp(fullData.lowSNRMask) = NaN;
    p1 = pcolor(fullData.mTime, range * 1e-3, transpose(rcsTmp) / lcMean * 1e6);
    p1.EdgeColor = 'none';
    caxis([0, 20]);
    colormap('jet');

    xlabel('时间');
    ylabel('距离 (千米)');
    title(sprintf('激光测雾雷达时空信号图 (%s)', datestr(fullData.mTime(1), 'yyyy-mm-dd')));

    xlim([fullData.mTime(1), fullData.mTime(end)]);
    ylim([0, 4]);

    set(gca, 'XMinorTick', 'on', 'xtick', linspace(fullData.mTime(1), fullData.mTime(end), 6), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 11);
    datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

    cb = colorbar('Position', [0.89, 0.2, 0.03, 0.6], 'Units', 'Normalized');
    titleHandle = get(cb, 'Title');
    set(titleHandle, 'String', '[Mm-1sr-1]');
    set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);
end

% full data ext
if exist('fullData', 'var')

    figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

    subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
    hold on;
    extTmp = fullData.ext;
    extTmp(fullData.lowSNRMask) = NaN;
    p1 = pcolor(fullData.mTime, range * 1e-3, transpose(extTmp) * 1e3);
    p1.EdgeColor = 'none';
    caxis([0, 10]);
    colormap('jet');

    xlabel('时间');
    ylabel('距离 (千米)');
    title(sprintf('消光系数 (%s)', datestr(fullData.mTime(1), 'yyyy-mm-dd')));

    xlim([fullData.mTime(1), fullData.mTime(end)]);
    ylim([0, 4]);

    set(gca, 'XMinorTick', 'on', 'xtick', linspace(fullData.mTime(1), fullData.mTime(end), 6), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 11);
    datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

    cb = colorbar('Position', [0.89, 0.2, 0.03, 0.6], 'Units', 'Normalized');
    titleHandle = get(cb, 'Title');
    set(titleHandle, 'String', '[km-1]');
    set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);
end

% full data vis
if exist('fullData', 'var')

    figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

    subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
    hold on;
    visTmp = fullData.vis;
    visTmp(fullData.lowSNRMask) = NaN;
    p1 = pcolor(fullData.mTime, range * 1e-3, transpose(visTmp) * 1e-3);
    p1.EdgeColor = 'none';
    caxis([0, 30]);
    load('vis_colormap.mat');
    colormap(double(visColorbar) / 255);

    xlabel('时间');
    ylabel('距离 (千米)');
    title(sprintf('能见度 (%s)', datestr(fullData.mTime(1), 'yyyy-mm-dd')));

    xlim([fullData.mTime(1), fullData.mTime(end)]);
    ylim([0, 4]);

    set(gca, 'XMinorTick', 'on', 'xtick', linspace(fullData.mTime(1), fullData.mTime(end), 6), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 11);
    datetick(gca, 'x', 'HH:MM', 'keepticks', 'keeplimits');

    cb = colorbar('Position', [0.89, 0.2, 0.03, 0.6], 'Units', 'Normalized');
    titleHandle = get(cb, 'Title');
    set(titleHandle, 'String', '[km]');
    set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);
end

% full data comparison
if exist('fullData', 'var')
    figure('Position', [0, 10, 550, 300], 'Units', 'Pixels', 'Color', 'w');

    hold on;
    s1 = plot(visData.mTime, visData.vis * 1e-3, 'color', 'k', 'Marker', '.', 'MarkerFaceColor', 'k', 'markeredgecolor', 'k', 'DisplayName', '能见度仪');
    s2 = plot(fullData.mTime, fullData.vis(:, 80) * 1e-3, 'color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b', 'markeredgecolor', 'b', 'DisplayName', sprintf('激光测雾雷达(%3d米)', floor(range(80))));
    hold off;

    xlim([fullData.mTime(1), fullData.mTime(end)]);
    ylim([0, 50]);

    xlabel('时间');
    ylabel('能见度 (千米)');
    title(sprintf('能见度序列对比(%s)', datestr(fullData.mTime(1), 'yyyy-mm-dd')));

    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);
    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');
    legend([s1, s2], 'location', 'NorthEast');
end