%------------------------------------------------------------------------------%
% Generate quicklooks from raw data
% 2023-05-05
% Zhenping Yin
% zp.yin@whu.edu.cn
%------------------------------------------------------------------------------%

%% Parameter Definition
dataPath = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\RawData';
savePath = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\Quicklooks';
location = 'QingDao1';
tRange = [datenum(2023, 12, 25, 0, 0, 0), datenum(2023, 12, 31, 23, 59, 59)];
visRetMethod = 'xian';   % xian: Xian's method; quasi: Quasi retrieval
debug = false;
olHeight = 1000;
overlapCor = true;
overlapFile = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\wxzk\near-sea-shore-measurement\2023-12-29\overlap_20231229.mat';

for iDate = floor(tRange(1)):floor(tRange(2))
    fprintf('Finished %6.2f%%\n', (iDate - floor(tRange(1))) / (floor(tRange(2)) - floor(tRange(1))) * 100);

    measDay = iDate;
    dayPath = fullfile(dataPath, ...
                       location, ...
                       sprintf('%d', str2double(datestr(measDay, 'yyyy'))), ...
                       sprintf('%d', str2double(datestr(measDay, 'mm'))), ...
                       sprintf('%d', str2double(datestr(measDay, 'dd'))));

    if ~ exist(dayPath, 'dir')
        continue;
    end

    %% Find All Data SubFolders
    allSubFolders = listdir(dayPath, '^H\w*', 2);

    %% Filter Data SubFolders
    subFolders = {};
    granuleTimes = [];
    for iGranule = 1:length(allSubFolders)
        tmp = basedir(allSubFolders{iGranule});
        thisGranuleTime = datenum(tmp((end - 13):end), 'yyyymmddHHMMSS');

        if (thisGranuleTime >= tRange(1)) && (thisGranuleTime <= tRange(2))
            subFolders = cat(2, subFolders, allSubFolders{iGranule});
            granuleTimes = cat(2, granuleTimes, thisGranuleTime);
        end
    end

    %% Iterate Over All subFolder
    im1 = [];
    im2 = [];
    for iFolder = 1:length(subFolders)
        fprintf('Finished %6.2f%%: processing %s\n', (iFolder - 1) / length(subFolders) * 100, subFolders{iFolder});

        %% Read Data
        dataFiles = listfile(subFolders{iFolder}, '\w*.VIS', 1);
        data = readVIS(dataFiles, 'debug', debug);

        %% Display RCS
        range = ((1:data.nBins(1)) + 0.5) * data.hRes(1) - 48.75;
        X = (sin(data.zenithAng' / 180 * pi) .* sin(data.azimuthAng' / 180 * pi)) * range;
        Y = (sin(data.zenithAng' / 180 * pi) .* cos(data.azimuthAng' / 180 * pi)) * range;
        bg = nanmean(squeeze(data.rawSignal(:, 1, 2900:2950)), 2);
        signal = squeeze(data.rawSignal(:, 1, :)) - repmat(bg, 1, data.nBins(1));
        if overlapCor
            ol = load(overlapFile);
            signal = signal ./ repmat(ol.ov', size(signal, 1), 1);
        end
        rcs = signal .* repmat(range, length(data.hRes), 1).^2;
        snr = (signal) ./ sqrt(squeeze(data.rawSignal(:, 1, :)));
        lowSNRMask = false(size(signal));
        for iPrf = 1:size(rcs, 1)
            snrPrf = snr(iPrf, :);
            snrPrf(range < olHeight) = NaN;
            snrLow = find(snrPrf < 1, 1);
            lowSNRMask(iPrf, snrLow:end) = true;
        end

        %% Visbility Retrieval
        ext = NaN(size(signal));
        for iPrf = 1:length(data.hRes)

            if strcmpi(visRetMethod, 'xian')
                ext(iPrf, :) = extRet_Xian(range, signal(iPrf, :), bg(iPrf), 'minSNR', 0.5, 'rangeFullOverlap', 200);
            elseif strcmpi(visRetMethod, 'quasi')
                [~, ext(iPrf, :)] = extRet_Holger(range, signal(iPrf, :), ...
                    'calibration_constant', 1.9e15, ...
                    'fullOverlapR', 280, ...
                    'elevation_angle', data.zenithAng(iPrf));
            else
            end

        end
        %ext(:, range <= 500) = NaN;

        %% Create Save Path
        subSavePath = fullfile(savePath, location, datestr(granuleTimes(iFolder), 'yyyy'), datestr(granuleTimes(iFolder), 'mm'), datestr(granuleTimes(iFolder), 'dd'));
        if ~ exist(subSavePath, 'dir')
            mkdir(subSavePath);
        end
        figure('Position', [0, 0, 700, 400], 'color', 'w', 'visible', 'off');
        rcs(lowSNRMask) = NaN;
        [~, p1] = polarPcolor(range / 1e3, data.azimuthAng, rcs, 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 10], 'Ncircles', 5, 'labelR', '', 'tickSize', 12, 'tickColor', 'm', 'typeRose', 'default', 'cRange', [0, 1e10]);
        ylabel(p1, '距离修正信号');
        set(p1, 'location', 'westoutside');
        colormap(gca, myColormap('jetImage'));
        text(0.3, 1.1, sprintf('%s', datestr(mean(data.startTime), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'Bold');
        text(0.6, -0.05, '距离 (千米)', 'Units', 'normalized', 'FontSize', 13, 'FontWeight', 'Bold');

        pngFile = fullfile(subSavePath, sprintf('%s_range_corrected_signal.png', datestr(granuleTimes(iFolder), 'yyyymmdd_HHMMSS')));
        saveas(gcf, pngFile, 'png');
        im1 = cat(4, im1, rgb2ind(imread(pngFile), myColormap('jetImage')));
        close;

       %% extinction
        % figure('Position', [0, 0, 700, 400], 'color', 'w', 'visible', 'off');
        % [~, p1] = polarPcolor(range / 1e3, data.azimuthAng, ext, 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 10], 'Ncircles', 5, 'labelR', '', 'tickSize', 12, 'tickColor', 'm', 'typeRose', 'default', 'cRange', [0, 1e-4]);
        % ylabel(p1, 'extinction (m-1)');
        % set(p1, 'location', 'westoutside');
        % colormap(gca, myColormap('jetImage'));
        % text(0.3, 1.1, sprintf('%s', datestr(mean(data.startTime), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'Bold');
        % text(0.6, -0.05, 'distance (km)', 'Units', 'normalized', 'FontSize', 13, 'FontWeight', 'Bold');

        % pngFile = fullfile(subSavePath, sprintf('%s_extinction_%s.png', datestr(granuleTimes(iFolder), 'yyyymmdd_HHMMSS'), visRetMethod));
        % saveas(gcf, pngFile, 'png');
        % im2 = cat(4, im2, rgb2ind(imread(pngFile), myColormap('jetImage')));
        % close;

        %% visiblity
        ext(lowSNRMask & (snr < 3)) = NaN;
        vis = ext2vis(ext);
        vis(isnan(vis)) = 1e5;

        figure('Position', [0, 0, 700, 400], 'color', 'w', 'visible', 'off');
        [~, p1] = polarPcolor(range / 1e3, data.azimuthAng, vis * 1e-3, 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 10], 'Ncircles', 5, 'labelR', '', 'tickSize', 12, 'tickColor', 'm', 'typeRose', 'default', 'cRange', [0, 30]);
        ylabel(p1, '能见度 (千米)');
        set(p1, 'location', 'westoutside');
        load('vis_colormap.mat');
        colormap(gca, double(visColorbar) / 255);
        text(0.3, 1.1, sprintf('%s', datestr(mean(data.startTime), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'Bold');
        text(0.6, -0.05, '距离 (千米)', 'Units', 'normalized', 'FontSize', 13, 'FontWeight', 'Bold');

        pngFile = fullfile(subSavePath, sprintf('%s_visibility_%s.png', datestr(granuleTimes(iFolder), 'yyyymmdd_HHMMSS'), visRetMethod));
        saveas(gcf, pngFile, 'png');
        im2 = cat(4, im2, rgb2ind(imread(pngFile), flipud(myColormap('jetImage'))));
        close;
    end

    imwrite(uint8(im1), myColormap('jetImage'), fullfile(subSavePath, sprintf('%s_rcs_animation.gif', datestr(mean(data.startTime), 'yyyymmdd'))), 'DelayTime', 0.3, 'LoopCount', inf, 'BackgroundColor', 1);
    imwrite(uint8(im2), flipud(myColormap('jetImage')), fullfile(subSavePath, sprintf('%s_vis_animation.gif', datestr(mean(data.startTime), 'yyyymmdd'))), 'DelayTime', 0.3, 'LoopCount', inf, 'BackgroundColor', 254);
end
