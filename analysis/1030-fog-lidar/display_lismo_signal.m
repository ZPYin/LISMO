clc;
global LISMO_VARS;

%% Parameter Initialization
dataFolder = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar';
date = datenum(2023, 6, 12);
hRes = 7.5;
firstRangeBin = 15;
deadtime = 22;   % [ns]

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

% deadtime correction
rawSignal = data.rawSignal / (50 * thisData.nShots(1) * 1e-3);
corSignal = rawSignal ./ (1 - deadtime * rawSignal * 1e-3);

bg = nanmean(corSignal(3500:3600, :, :), 1);
signal = corSignal - repmat(bg, size(corSignal, 1), 1, 1);
height = ((1:size(corSignal, 1)) - firstRangeBin + 0.5) * hRes;
rcs = signal .* repmat(reshape(height, length(height), 1, 1), 1, size(signal, 2), size(signal, 3)).^2;

%% Data Display
iProf = 665;

% profile with full range
figure('Position', [0, 0, 500, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(1:size(corSignal, 1), squeeze(corSignal(:, 2, iProf)), '-', 'Color', [15, 15, 172] / 255, 'LineWidth', 2, 'DisplayName', 'Far-Range');
p2 = plot(1:size(corSignal, 1), squeeze(corSignal(:, 1, iProf)), '-', 'Color', [56, 87, 35] / 255, 'LineWidth', 2, 'DisplayName', 'Near-Range');

xlim([1, size(corSignal, 1)]);

xlabel('Range Bin');
ylabel('Raw Signal (MHz)');

title(sprintf('Time: %s', datestr(data.mTime(iProf), 'yyyy-mm-dd HH:MM:SS')));

set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log');

legend([p1, p2], 'Location', 'NorthEast');

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'test-measurement-single-profile.png'), '-r300');

% profile with limited range
figure('Position', [0, 0, 500, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(1:size(corSignal, 1), squeeze(corSignal(:, 2, iProf)), '-', 'Color', [15, 15, 172] / 255, 'LineWidth', 2, 'DisplayName', 'Far-Range');
p2 = plot(1:size(corSignal, 1), squeeze(corSignal(:, 1, iProf)), '-', 'Color', [56, 87, 35] / 255, 'LineWidth', 2, 'DisplayName', 'Near-Range');

xlim([1, 200]);

xlabel('Range Bin');
ylabel('Raw Signal (MHz)');

title(sprintf('Time: %s', datestr(data.mTime(iProf), 'yyyy-mm-dd HH:MM:SS')));

set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log');

legend([p1, p2], 'Location', 'NorthEast');

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'test-measurement-single-profile-zoomin.png'), '-r300');

% range corrected signal
figure('Position', [0, 0, 500, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(height, nanmean(squeeze(corSignal(:, 2, iProf)), 2), '-', 'Color', [15, 15, 172] / 255, 'LineWidth', 2, 'DisplayName', 'Far-Range');
p2 = plot(height, nanmean(squeeze(corSignal(:, 1, iProf)), 2), '-', 'Color', [56, 87, 35] / 255, 'LineWidth', 2, 'DisplayName', 'Near-Range');

xlim([0, max(height)]);

xlabel('Distance [m]');
ylabel('range cor. sig. [a.u.]');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'linear', 'TickLength', [0.03, 0.03]);

legend([p1, p2], 'Location', 'NorthEast');

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'test-measurement-single-profile-rcs.png'), '-r300');

%% range corrected signal (near-range)
figure('Position', [0, 0, 500, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(height, nanmean(squeeze(rcs(:, 2, (iProf - 60):iProf)), 2), '-', 'Color', [15, 15, 172] / 255, 'LineWidth', 2, 'DisplayName', 'Far-Range');
p2 = plot(height, nanmean(squeeze(rcs(:, 1, (iProf - 60):iProf)), 2), '--', 'Color', [56, 87, 35] / 255, 'LineWidth', 1, 'DisplayName', 'Near-Range');
p3 = plot(height, nanmean(squeeze(rcs(:, 1, (iProf - 60):iProf)), 2) * 140000.0, '-', 'Color', [56, 87, 35] / 255, 'Marker', 'o', 'MarkerSize', 5, 'LineWidth', 2, 'DisplayName', 'Near-Range (x k)');

xlim([0, 500]);

xlabel('Distance [m]');
ylabel('range cor. sig. [a.u.]');

title(sprintf('Time: %s', datestr(data.mTime(iProf), 'yyyy-mm-dd HH:MM:SS')));

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log', 'TickLength', [0.03, 0.03]);

legend([p1, p2, p3], 'Location', 'SouthEast');

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'test-measurement-single-profile-rcs-near-range.png'), '-r300');