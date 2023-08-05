clc;
close all;

%% Parameter Definition
path1030 = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\1030vs1064\1030';
path1064 = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\1030vs1064\1064';
firstRangeBin1030 = 15;
firstRangeBin1064 = 15;
maxBin1030 = 3400;
hRes1030 = 7.5;
hRes1064 = 15;
deadtime = 0;
flagRead1030 = true;
flagRead1064 = true;
flagShowTHPlot = true;

%% Read Data

% 1030
if flagRead1030
    dataFiles = listfile(path1030, '.*.dat', 1);
    data1030 = struct();
    data1030.mTime = [];
    data1030.rawSignal = [];
    for iFile = 1:length(dataFiles)
        fprintf('Finished %6.2f%%: reading %s\n', (iFile - 1) / length(dataFiles) * 100, dataFiles{iFile});
        thisData = readALADat(dataFiles{iFile});

        data1030.mTime = cat(2, data1030.mTime, thisData.mTime);
        data1030.rawSignal = cat(3, data1030.rawSignal, reshape(thisData.rawSignal(1:maxBin1030, :), maxBin1030, size(thisData.rawSignal, 2), 1));
    end

    % deadtime correction
    data1030.rawSignal = data1030.rawSignal / (50 * thisData.nShots(1) * 1e-3);
    data1030.corSignal = data1030.rawSignal ./ (1 - deadtime * data1030.rawSignal * 1e-3);

    data1030.bg = nanmean(data1030.corSignal(3200:3300, :, :), 1);
    data1030.signal = data1030.corSignal - repmat(data1030.bg, size(data1030.corSignal, 1), 1, 1);
    data1030.height = ((1:size(data1030.corSignal, 1)) - firstRangeBin1030 + 0.5) * hRes1030;
    data1030.rcs = data1030.signal .* repmat(reshape(data1030.height, length(data1030.height), 1, 1), 1, size(data1030.signal, 2), size(data1030.signal, 3)).^2;
end

% 1064
if flagRead1064
    if ~ isempty(path1064)
        dataFiles = listfile(path1064, '.*.dat', 1);
    else
        dataFiles = {};
    end
    data1064 = struct();
    data1064.mTime = [];
    data1064.rawSignal = [];
    for iFile = 1:length(dataFiles)
        fprintf('Finished %6.2f%%: reading %s\n', (iFile - 1) / length(dataFiles) * 100, dataFiles{iFile});
        thisData = readALADat(dataFiles{iFile});

        data1064.mTime = cat(2, data1064.mTime, thisData.mTime);
        data1064.rawSignal = cat(3, data1064.rawSignal, reshape(thisData.rawSignal, size(thisData.rawSignal, 1), size(thisData.rawSignal, 2), 1));
    end

    % deadtime correction
    if ~ isempty(dataFiles)
        data1064.rawSignal = data1064.rawSignal / (50 * thisData.nShots(1) * 1e-3);
        data1064.corSignal = data1064.rawSignal ./ (1 - deadtime * data1064.rawSignal * 1e-3);

        data1064.bg = nanmean(data1064.corSignal(3200:3300, :, :), 1);
        data1064.signal = data1064.corSignal - repmat(data1064.bg, size(data1064.corSignal, 1), 1, 1);
        data1064.height = ((1:size(data1064.corSignal, 1)) - firstRangeBin1064 + 0.5) * hRes1064;
        data1064.rcs = data1064.signal .* repmat(reshape(data1064.height, length(data1064.height), 1, 1), 1, size(data1064.signal, 2), size(data1064.signal, 3)).^2;
    else
        data1064.rawSignal = [];
        data1064.corSignal = [];
        data1064.bg = [];
        data1064.signal = [];
        data1064.height = [];
        data1064.rcs = [];
    end
end

%% Data Visualization

if flagShowTHPlot
    figure('Position', [0, 0, 500, 500], 'Units', 'Pixels', 'Color', 'w');

    subplot(211);
    rcs = squeeze(data1030.rcs(:, 2, :));
    rcs(rcs <= 0) = NaN;
    p1 = pcolor(data1030.mTime, data1030.height / 1e3, log10(rcs));
    p1.EdgeColor = 'none';
    caxis([3, 6]);
    colormap('jet');

    xlabel('');
    ylabel('Distance (km)');
    title('1030 nm lidar (elevation angle 7\circ)');
    
    xlim([datenum(2023, 6, 26, 19, 0, 0), datenum(2023, 6, 26, 22, 30, 0)]);
    ylim([0, 15]);

    colorbar();

    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Tickdir', 'out', 'Layer', 'top', 'box', 'on');
    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    subplot(212);
    rcs = squeeze(data1064.rcs(:, 2, :));
    rcs(rcs <= 0) = NaN;
    p1 = pcolor(data1064.mTime, data1064.height / 1e3, log10(rcs));
    p1.EdgeColor = 'none';
    caxis([3, 6.5]);
    colormap('jet');

    xlabel('');
    ylabel('Distance (km)');
    title('1064 nm lidar (elevation angle 7\circ)');
    
    xlim([datenum(2023, 6, 26, 19, 0, 0), datenum(2023, 6, 26, 22, 30, 0)]);
    ylim([0, 15]);

    colorbar();

    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Tickdir', 'out', 'Layer', 'top', 'box', 'on');
    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');
end

% focus on blind zone
tRangeBZ = [datenum(2023, 6, 26, 22, 12, 0), datenum(2023, 6, 26, 22, 22, 0)];
isInTRange1030 = (data1030.mTime >= tRangeBZ(1)) & (data1030.mTime <= tRangeBZ(2));
isInTRange1064 = (data1064.mTime >= tRangeBZ(1)) & (data1064.mTime <= tRangeBZ(2));

% range corrected signal
figure('Position', [0, 0, 500, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(data1030.height, nanmean(squeeze(data1030.rcs(:, 2, isInTRange1030)), 2), '-', 'Color', [15, 15, 172] / 255, 'LineWidth', 2, 'DisplayName', '1030 FR');
p2 = plot(data1030.height, nanmean(squeeze(data1030.rcs(:, 1, isInTRange1030)), 2), '-', 'Color', [56, 87, 35] / 255, 'LineWidth', 2, 'DisplayName', '1030 NR');
if ~ isempty(data1064.corSignal)
    p3 = plot(data1064.height, nanmean(squeeze(data1064.rcs(:, 1, isInTRange1064)), 2), '--', 'Color', 'r', 'LineWidth', 2, 'DisplayName', '1064');
else
    p3 = plot([], [], '--', 'Color', 'r', 'LineWidth', 2, 'DisplayName', '1064');
end
hold off;

xlim([0, max(data1030.height)]);

xlabel('Distance [m]');
ylabel('range cor. sig. [a.u.]');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'linear', 'TickLength', [0.03, 0.03]);

legend([p1, p2, p3], 'Location', 'NorthEast');

% export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'test-measurement-single-profile-rcs.png'), '-r300');

%% range corrected signal (near-range)
figure('Position', [0, 0, 500, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(data1030.height, nanmean(squeeze(data1030.rcs(:, 2, isInTRange1030)), 2), '-', 'Color', [15, 15, 172] / 255, 'LineWidth', 2, 'DisplayName', '1030 FR');
p2 = plot(data1030.height, nanmean(squeeze(data1030.rcs(:, 1, isInTRange1030)), 2), '--', 'Color', [56, 87, 35] / 255, 'LineWidth', 1, 'DisplayName', '1030 NR');
p3 = plot(data1030.height, nanmean(squeeze(data1030.rcs(:, 1, isInTRange1030)), 2) * 140000.0, '-', 'Color', [56, 87, 35] / 255, 'Marker', 'o', 'MarkerSize', 5, 'LineWidth', 2, 'DisplayName', 'Near-Range (x k)');
if ~ isempty(data1064.rcs)
    p4 = plot(data1064.height, nanmean(squeeze(data1064.rcs(:, 1, isInTRange1064)), 2), '--', 'Color', 'r', 'LineWidth', 2, 'DisplayName', '1064');
else
    p4 = plot([], [], '--', 'Color', 'r', 'LineWidth', 2, 'DisplayName', '1064');
end
hold off;

xlim([0, 700]);

xlabel('Distance [m]');
ylabel('range cor. sig. [a.u.]');

title(sprintf('Time: %s', datestr(mean(data1030.mTime(isInTRange1030)), 'yyyy-mm-dd HH:MM:SS')));

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YScale', 'log', 'TickLength', [0.03, 0.03]);

legend([p1, p2, p3, p4], 'Location', 'SouthEast');