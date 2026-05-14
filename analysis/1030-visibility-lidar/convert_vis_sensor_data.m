% convert visibility sensor data.

%% Parameter Definition
visFile = 'G:\backup\vis-lidar\20260428-tianjin-evaluation\前散仪对比数据\2026年05月01日00时+至+2026年05月07日00时+观测试验基地+____+能见度对比+分钟数据.xls';

%% Read Data
[~, ~, dataCell] = xlsread(visFile);

%% Data Conversion
nRec = size(dataCell, 1) - 1;
mTime = NaN(1, nRec);
vis = NaN(1, nRec);
for iR = 1:nRec
    mTime(iR) = datenum(dataCell{iR + 1, 2}, 'yyyy/mm/dd HH:MM:SS');
    vis(iR) = dataCell{iR + 1, 4};
end
visData.mTime = mTime;
visData.vis = vis;

%% Save Data
save('vis-sensor-data.mat', 'vis');

%% Display
figure('Position', [0, 30, 500, 250], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(visData.mTime, visData.vis * 1e-3, '-r', 'DisplayName', dataCell{2, 1});
hold off;

xlim([min(visData.mTime), max(visData.mTime)]);
ylim([0, 50]);

xlabel('时间');
ylabel('能见度 (千米)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);
datetick(gca, 'x', 'mm-dd', 'keeplimits', 'keepticks');

legend([p1], 'Location', 'NorthWest');

export_fig(gcf, 'vis_sensor_overview.png', '-r300');