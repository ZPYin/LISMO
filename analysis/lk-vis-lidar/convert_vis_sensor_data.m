% convert visibility sensor data.

%% Parameter Definition
visFile = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\能见度仪分钟数据_1101--1128.xlsx';

%% Read Data
[~, ~, data_998] = xlsread(visFile, 'I0998.txt');
[~, ~, data_297] = xlsread(visFile, 'I0297.txt');

%% Data Conversion
nRecord1 = size(data_998, 1) - 1;
mTime1 = NaN(1, nRecord1);
vis1 = NaN(1, nRecord1);
for iR = 1:nRecord1
    mTime1(iR) = datenum(data_998{iR + 1, 1}, 'yyyy-mm-dd HH:MM:SS.FFF');
    vis1(iR) = str2double(data_998{iR + 1, 2});
end
vis998.mTime = mTime1;
vis998.vis = vis1;

nRecord2 = size(data_297, 1) - 1;
mTime2 = NaN(1, nRecord2);
vis2 = NaN(1, nRecord2);
for iR = 1:nRecord2
    mTime2(iR) = datenum(data_297{iR + 1, 1}, 'yyyy-mm-dd HH:MM:SS.FFF');
    vis2(iR) = str2double(data_297{iR + 1, 2});
end
vis297.mTime = mTime2;
vis297.vis = vis2;

%% Save Data
save('vis-sensor-data.mat', 'vis998', 'vis297');

%% Display
figure('Position', [0, 30, 500, 250], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(vis297.mTime, vis297.vis * 1e-3, '-r', 'DisplayName', 'I0297');
p2 = plot(vis998.mTime, vis998.vis * 1e-3, '-k', 'DisplayName', 'I0998');
hold off;

xlim([min(vis297.mTime), max(vis297.mTime)]);
ylim([0, 50]);

xlabel('时间');
ylabel('能见度 (千米)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);
datetick(gca, 'x', 'mm-dd', 'keeplimits', 'keepticks');

legend([p1, p2], 'Location', 'NorthWest');

export_fig(gcf, 'vis_sensor_overview.png', '-r300');