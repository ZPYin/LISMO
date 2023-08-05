clc;

%% Parameter Initialization
dataFolder = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\data\1030';
saveFolder = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\quicklooks';
date = datenum(2023, 7, 25);
hRes = 7.5;
firstRangeBin = 15;
cRange = [0, 1e6];
flagReadData = true;
deadtime = 22;   % [ns]

%% List Data Files
dataFiles = listfile(fullfile(dataFolder, datestr(date, 'yyyy-mm-dd')), '.*dat', 1);

%% Read Lidar Data
if flagReadData
    data = struct;
    data.mTime = [];
    data.rawSignal = [];
    for iFile = 1:length(dataFiles)
        fprintf('Finished %6.2f%%: reading %s\n', (iFile - 1) / length(dataFiles) * 100, dataFiles{iFile});
        thisData = readALADat(dataFiles{iFile});
        
        data.mTime = cat(2, data.mTime, thisData.mTime);
        data.rawSignal = cat(3, data.rawSignal, reshape(thisData.rawSignal, size(thisData.rawSignal, 1), size(thisData.rawSignal, 2), 1));
    end

    % deadtime correction
    rawSignal = data.rawSignal / (50 * thisData.nShots(1) * 1e-3);
    corSignal = rawSignal ./ (1 - deadtime * rawSignal * 1e-3);

    bg = nanmean(corSignal(3000:3100, :, :), 1);
    signal = corSignal - repmat(bg, size(corSignal, 1), 1, 1);
    height = ((1:size(corSignal, 1)) - firstRangeBin + 0.5) * hRes;
    rcs = signal .* repmat(reshape(height, length(height), 1, 1), 1, size(signal, 2), size(signal, 3)).^2;
end

%% Data Display
figure('Position', [0, 0, 500, 250], 'Units', 'Pixels', 'Color', 'w');

p1 = pcolor(data.mTime, height, squeeze(rcs(:, 2, :)));
p1.EdgeColor = 'None';

xlabel('Time (HH:MM)');
ylabel('Distance (m)');
title(sprintf('Range corrected signal at %s', datestr(data.mTime(1), 'yyyy-mm-dd')));

caxis(cRange);
colormap('jet');

set(gca, 'TickDir', 'out', 'Box', 'on');
xlim([min(data.mTime), max(data.mTime)]);
ylim([0, 20000]);
datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

colorbar();

export_fig(gcf, fullfile(saveFolder, sprintf('%s_rcs.png', datestr(data.mTime(1), 'yyyy-mm-dd'))), '-r300');