clc;

%% Parameter Initialization
dataFolder = 'C:\Users\ZPYin\Documents\Data\Airda-Data\2023\Q2';
date = datenum(2023, 6, 13);
hRes = 15;
firstRangeBin = 10;
deadtime = 1;   % [ns]
iProf = 480:550;

%% List Data Files
dataFiles = listfile(fullfile(dataFolder, datestr(date, 'yyyy-mm-dd')), '.*dat', 1);

%% Read Lidar Data
data = struct;
data.mTime = [];
data.rawSignal = [];
for iFile = 1:length(dataFiles)
    fprintf('Finished %6.2f%%: reading %s\n', (iFile - 1) / length(dataFiles) * 100, dataFiles{iFile});
    thisData = readALADat(dataFiles{iFile});

    data.mTime = cat(2, data.mTime, thisData.mTime);
    data.rawSignal = cat(3, data.rawSignal, reshape(thisData.rawSignal, size(thisData.rawSignal, 1), size(thisData.rawSignal, 2), 1));
end

%% deadtime correction
rawSignal = data.rawSignal / (50 * thisData.nShots(1) * 1e-3);
corSignal = rawSignal ./ (1 - deadtime * rawSignal * 1e-3);

bg = nanmean(corSignal(1950:2000, :, :), 1);
signal = corSignal - repmat(bg, size(corSignal, 1), 1, 1);
height = ((1:size(corSignal, 1)) - firstRangeBin + 0.5) * hRes;
rcs = signal .* repmat(reshape(height, length(height), 1, 1), 1, size(signal, 2), size(signal, 3)).^2;

%% Visualization
figure('Position', [0, 0, 600, 500], 'Units', 'Pixels', 'Color', 'w');

subplot(211);
lineInstances = [];
hold on;
for iCh = 1:size(signal, 2)
    p1 = plot(height / 1e3, squeeze(nanmean(signal(:, iCh, iProf), 3)), 'DisplayName', thisData.channelLabel{iCh});
    lineInstances = cat(2, lineInstances, p1);
end

xlabel('Height (km)');
ylabel('Signal (MHz)');
title(sprintf('%s', datestr(date, 'yyyy-mm-dd')));

xlim([0, 10]);
ylim([1e-5, 1e2]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log');
legend(lineInstances, 'Location', 'NorthEast');

subplot(212);
lineInstances = [];
hold on;
for iCh = 1:size(signal, 2)
    p1 = plot(height / 1e3, squeeze(nanmean(rcs(:, iCh, iProf), 3)), 'DisplayName', thisData.channelLabel{iCh});
    lineInstances = cat(2, lineInstances, p1);
end

xlabel('Height (km)');
ylabel('RCS (a.u.)');
title(sprintf('%s', datestr(date, 'yyyy-mm-dd')));

xlim([0, 10]);
ylim([1e2, 1e8]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log');
legend(lineInstances, 'Location', 'NorthEast');

