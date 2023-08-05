%------------------------------------------------------------------------------%
% Generate quicklooks from raw data
% 2023-05-05
% Zhenping Yin
% zp.yin@whu.edu.cn
%------------------------------------------------------------------------------%

%% Parameter Definition
dataPath = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\RawData';
savePath = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\Quicklooks';
location = 'QingDao';
tRange = [datenum(2023, 7, 3, 0, 0, 0), datenum(2023, 7, 3, 23, 59, 59)];
visRetMethod = 'quasi';   % xian: Xian's method; quasi: Quasi retrieval
debug = false;
overlapCor = false;

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
            load('overlap.mat');
            signal = signal ./ repmat(ov, size(signal, 1), 1);
        end
        rcs = signal .* repmat(range, length(data.hRes), 1).^2;
        snr = (signal) ./ sqrt(squeeze(data.rawSignal(:, 1, :)));

        %% Visbility Retrieval
        ext = NaN(size(signal));
        for iPrf = 1:length(data.hRes)

            if strcmpi(visRetMethod, 'xian')
                ext(iPrf, :) = extRet_Xian(range, signal(iPrf, :), bg(iPrf), 'minSNR', 0.1);
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
        
        figure('color', 'w', 'visible', 'off');
        rcs(snr <= 1) = NaN;
        [~, p1] = polarPcolor(range / 1e3, data.azimuthAng, rcs, 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 10], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 4e9]);
        text(0.3, 1.2, sprintf('%s', datestr(mean(data.startTime), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'Bold');
        ylabel(p1, 'range cor. sig.');
        set(p1, 'location', 'westoutside');
        colormap(gca, myColormap('jetImage'));
        text(0.6, -0.15, 'distance (km)', 'Units', 'normalized', 'FontSize', 11, 'FontWeight', 'light');

        pngFile = fullfile(subSavePath, sprintf('%s_range_corrected_signal.png', datestr(granuleTimes(iFolder), 'yyyymmdd_HHMMSS')));
        export_fig(gcf, pngFile, '-r300');
        im1 = cat(4, im1, rgb2ind(imread(pngFile), myColormap('jetImage')));
        close;

       %% extinction
        figure('color', 'w', 'visible', 'on');
        [~, p1] = polarPcolor(range / 1e3, data.azimuthAng, ext, 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 10], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 1e-4]);
        text(0.3, 1.2, sprintf('%s', datestr(mean(data.startTime), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'Bold');
        ylabel(p1, 'extinction (m-1)');
        set(p1, 'location', 'westoutside');
        colormap(gca, myColormap('jetImage'));
        text(0.6, -0.15, 'distance (km)', 'Units', 'normalized', 'FontSize', 11, 'FontWeight', 'light');

        pngFile = fullfile(subSavePath, sprintf('%s_extinction_%s.png', datestr(granuleTimes(iFolder), 'yyyymmdd_HHMMSS'), visRetMethod));
        export_fig(gcf, pngFile, '-r300');
        im2 = cat(4, im2, rgb2ind(imread(pngFile), myColormap('jetImage')));
        close;

        %% Extinction Retrieval
        ext(snr < 3) = NaN;
        vis = ext2vis(ext);
        vis(isnan(vis)) = 1e5;

        figure('color', 'w', 'visible', 'on');
        [~, p1] = polarPcolor(range / 1e3, data.azimuthAng, vis, 'Nspokes', 7, 'colormap', 'hot', 'GridLineStyle', '--', 'RLim', [0, 10], 'Ncircles', 5, 'labelR', '', 'typeRose', 'default', 'cRange', [0, 2e4]);
        text(0.3, 1.2, sprintf('%s', datestr(mean(data.startTime), 'yyyy-mm-dd HH:MM')), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'Bold');
        ylabel(p1, 'visiblity (m)');
        set(p1, 'location', 'westoutside');
        colormap(gca, flipud(myColormap('jetImage')));
        text(0.6, -0.15, 'distance (km)', 'Units', 'normalized', 'FontSize', 11, 'FontWeight', 'light');

        pngFile = fullfile(subSavePath, sprintf('%s_visibility_%s.png', datestr(granuleTimes(iFolder), 'yyyymmdd_HHMMSS'), visRetMethod));
        export_fig(gcf, pngFile, '-r300');
        im2 = cat(4, im2, rgb2ind(imread(pngFile), flipud(myColormap('jetImage'))));
        close;
    end

    imwrite(uint8(im1), myColormap('jetImage'), fullfile(subSavePath, sprintf('%s_rcs_animation.gif', datestr(mean(data.startTime), 'yyyymmdd'))), 'DelayTime', 0.3, 'LoopCount', inf, 'BackgroundColor', 1);
    imwrite(uint8(im2), flipud(myColormap('jetImage')), fullfile(subSavePath, sprintf('%s_vis_animation.gif', datestr(mean(data.startTime), 'yyyymmdd'))), 'DelayTime', 0.3, 'LoopCount', inf, 'BackgroundColor', 254);
end
