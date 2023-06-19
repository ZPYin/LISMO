clc;
addpath('C:\Users\ZPYin\Documents\Coding\Matlab\Radiosonde\matlab');

%% Parameter Initialization
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_mwl\设备更改后数据\03 1期 甘肃酒泉三波长\甘肃酒泉6.18数据';
% date = datenum(2023, 6, 13);
hRes = 15;
firstRangeBin = 17;
deadtime = 1;   % [ns]
flagReadData = true;
flagReadRS = true;

%% List Data Files
dataFiles = listfile(fullfile(dataFolder), '.*dat', 1);

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
end

%% deadtime correction
rawSignal = data.rawSignal / (thisData.nShots(1)) * 300 / thisData.hRes(1);
corSignal = rawSignal ./ (1 - deadtime * rawSignal * 1e-3);

bg = nanmean(corSignal(1950:2000, :, :), 1);
signal = corSignal - repmat(bg, size(corSignal, 1), 1, 1);
height = ((1:size(corSignal, 1)) - firstRangeBin + 0.5) * hRes;
rcs = signal .* repmat(reshape(height, length(height), 1, 1), 1, size(signal, 2), size(signal, 3)).^2;

%% Read Radiosonde Data
[alt, temp, pres, rh] = read_websonde(datenum(2023, 6, 18, 12, 0, 0), [datenum(2023, 6, 18, 11, 0, 0), datenum(2023, 6, 18, 14, 0, 0)], 52533, 'BUFR');
es = saturated_vapor_pres(temp);
WVMR_rs = 18.0160 / 28.9660 * rh / 100 .* es ./ (pres - rh / 100 .* es) * 1000;

%% Visualization
iProf = 1:30;
figure('Position', [0, 0, 600, 500], 'Units', 'Pixels', 'Color', 'w');

subplot(211);
lineInstances = [];
hold on;
for iCh = 1:size(signal, 2)
    p1 = plot(height / 1e3, squeeze(nanmean(signal(:, iCh, iProf), 3)), 'LineWidth', 2, 'DisplayName', thisData.channelLabel{iCh});
    lineInstances = cat(2, lineInstances, p1);
end

xlabel('Height (km)');
ylabel('Signal (MHz)');
title(sprintf('%s', datestr(data.mTime(1), 'yyyy-mm-dd HH:MM')));

xlim([0, 10]);
ylim([1e-5, 1e2]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log');
legend(lineInstances, 'Location', 'NorthEast');

subplot(212);
lineInstances = [];
hold on;
for iCh = 1:size(signal, 2)
    p1 = plot(height / 1e3, squeeze(nanmean(rcs(:, iCh, iProf), 3)), 'LineWidth', 2, 'DisplayName', thisData.channelLabel{iCh});
    lineInstances = cat(2, lineInstances, p1);
end

xlabel('Height (km)');
ylabel('RCS (a.u.)');
title(sprintf('%s', datestr(data.mTime(1), 'yyyy-mm-dd HH:MM')));

xlim([0, 10]);
ylim([1e2, 1e8]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log');
legend(lineInstances, 'Location', 'NorthEast');

figure('Position', [0, 0, 400, 500], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(mean(signal(:, 4, iProf), 3) ./ mean(signal(:, 3, iProf), 3) * 75 * 3.4 / 5, height, '-b', 'LineWidth', 2, 'DisplayName', 'lidar');
p2 = plot(WVMR_rs, alt - alt(1), '--r', 'LineWidth', 2, 'DisplayName', 'radiosonde');

xlabel('水汽混合比 (g/kg)');
ylabel('高度 (m)');

xlim([0, 15]);
ylim([0, 10000]);

legend([p1, p2], 'Location', 'NorthEast');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');