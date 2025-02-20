% 反演得到中科光博能见度激光雷达反演结果
% 作者：殷振平
% 日期：2025年2月9日

clc;
close all;

%% Parameter Definition
dataPath = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\ZKGB\data\2025-01-01';
resPath = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\ZKGB\results\L1';
figPath = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\ZKGB\figures';
flagOLCor = false;
olFile = 'overlap_20250101.mat';
AEConvFactor = (1064/550) ^1;
hFullOL = 500;
distOffset = -15;
visible = 'off';

%% Find Scans
scanPaths = listdir(dataPath, '\w{14}', 1);
for iScan = 1:length(scanPaths)
    fprintf('Finished %6.2f%%: reading %s\n', (iScan - 1) / length(scanPaths) * 100, scanPaths{iScan});

    thisScan = scanPaths{iScan};

    %% Read Data

    thisData = readVIS_ZKGB(fullfile(thisScan, sprintf('%s_1064Signal.xlsx', basename(thisScan))), 'debug', false);

    %% Preprocessing
    range = thisData.range + distOffset;
    bg = mean(thisData.signal(:, (end - 20):end), 2);
    noise = std(thisData.signal(:, (end - 20):end), 0, 2);
    signal = thisData.signal - repmat(bg, 1, thisData.nBins);
    angle = thisData.angle;
    mTime = thisData.mTime;
    if flagOLCor
        ol = load(olFile);
        signal = signal ./ repmat(transpose(ol.ov), size(thisData.signal, 1), 1);
    end
    rcs = signal .* repmat(range, size(thisData.signal, 1), 1).^2;
    snr = signal ./ repmat(noise, 1, thisData.nBins);

    %% Rayleigh Scattering
    [mBsc, mExt] = MolModel(ones(size(range)) * 0.005 * 1e3, 1064, 'meteor', 'standard_atmosphere');
    [~, mExt550] = MolModel(ones(size(range)) * 0.005 * 1e3, 550, 'meteor', 'standard_atmosphere');

    %% extinction retrieval
    extMat_Fernald = NaN(size(signal));
    for iPrf = 1:size(signal, 1)
        extMat_Fernald(iPrf, :) = extRet_Fernald(range, signal(iPrf, :), bg(iPrf), mBsc, mExt, 'snr', snr(iPrf, :), 'minSNR', 3, 'hFullOL', 500);
    end

    %% Extinction to Visibility
    visMat_Fernald = ext2vis(extMat_Fernald * AEConvFactor + mExt550);

    %% Save Results
    save(fullfile(resPath, sprintf('%s_L1.mat', basename(thisScan))), 'range', 'angle', 'mTime', 'extMat_Fernald', 'visMat_Fernald');

    %% Display

    % range corrected signal
    figure('Position', [100, 100, 600, 400], 'color', 'w', 'visible', visible);

    subplot('Position', [0.1, 0.1, 0.84, 0.84], 'Units', 'normalized');

    rcs(snr <= 3) = NaN;
    [~, p1] = polarPcolor(range / 1e3, transpose(angle), transpose(rcs), 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 5], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 0.5e7], 'tickSize', 12, 'tickColor', 'm');
    ylabel(p1, 'range cor. sig.');
    set(p1, 'location', 'westoutside', 'FontSize', 12);
    colormap(gca, myColormap('jetImage'));
    text(0.3, 1.2, sprintf('%s', datestr(mTime(1), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'Bold');
    text(0.6, -0.05, 'distance (km)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'Bold');

    % extinction
    figure('Position', [100, 100, 600, 400], 'color', 'w', 'visible', visible);

    subplot('Position', [0.1, 0.1, 0.84, 0.84], 'Units', 'normalized');

    ext_Fernald(snr <= 3) = NaN;
    [~, p1] = polarPcolor(range / 1e3, transpose(angle), transpose(extMat_Fernald) * 1e3, 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 5], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 0.5], 'tickSize', 12, 'tickColor', 'm');
    ylabel(p1, 'extinction (km-1)');
    set(p1, 'location', 'westoutside', 'FontSize', 12);
    colormap(gca, myColormap('jetImage'));
    text(0.3, 1.2, sprintf('%s', datestr(mTime(1), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'Bold');
    text(0.6, -0.05, 'distance (km)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'Bold');

    % visibility
    figure('Position', [100, 100, 600, 400], 'color', 'w', 'visible', visible);

    subplot('Position', [0.1, 0.1, 0.84, 0.84], 'Units', 'normalized');

    visMat_Fernald(snr <= 3) = NaN;
    [~, p1] = polarPcolor(range / 1e3, transpose(angle), transpose(visMat_Fernald) * 1e-3, 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 5], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 30], 'tickSize', 12, 'tickColor', 'm');
    ylabel(p1, 'visibility (km)');
    set(p1, 'location', 'westoutside', 'FontSize', 12);
    load('vis_colormap.mat');
    colormap(gca, double(visColorbar) / 255);
    text(0.3, 1.2, sprintf('%s', datestr(mTime(1), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'Bold');
    text(0.6, -0.05, 'distance (km)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'Bold');

    close all;

end