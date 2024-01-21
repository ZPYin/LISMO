clc;
close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\RawData\QingDao1\2023\12\26\H_SCAN15_135_2_20231226001201';
tRange = [datenum(2023, 12, 26, 0, 12, 26), datenum(2023, 12, 26, 0, 20, 0)];   % temporal range for data input
hRangeDisplay = [0, 200];
flagReadData = true;

%% Read Data
if flagReadData
    lData = readVIS(dataFolder, 'tRange', tRange, 'isDir', true, 'debug', true);
end

%% Preprocess
range = ((1:data.nBins(1)) - 0.5) * data.hRes(1) ;

%% Display
figure('Position', [0, 10, 500, 300], 'Units', 'Pixels', 'Color', 'w');

p1 = plot(range, squeeze(nanmean(lData.rawSignal(:, 1, :))), '-k', 'DisplayName', 'Ch 1'); hold on;
p2 = plot(range, squeeze(nanmean(lData.rawSignal(:, 2, :))), '-r', 'DisplayName', 'Ch 2'); hold on;

xlabel('Range (m)');
ylabel('Raw Signal');

xlim(hRangeDisplay);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);
legend([p1, p2], 'Location', 'NorthEast');