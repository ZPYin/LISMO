file = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\quicklooks\vis_statistics.txt';

mTime = [];
nProfiles = [];
meanErrFernald0297 = [];
stdErrFernald0297 = [];
outliersFernald0297 = [];
meanErrLK0297 = [];
stdErrLK0297 = [];
outliersLK0297 = [];
meanErrXian0297 = [];
stdErrXian0297 = [];
outliersXian0297 = [];

meanErrFernald0998 = [];
stdErrFernald0998 = [];
outliersFernald0998 = [];
meanErrLK0998 = [];
stdErrLK0998 = [];
outliersLK0998 = [];
meanErrXian0998 = [];
stdErrXian0998 = [];
outliersXian0998 = [];

fid = fopen(file, 'r');

while ~feof(fid)
    thisLine = fgetl(fid);
    if isempty(thisLine)
        break;
    end

    subStrs = strsplit(thisLine, ' ');
    mTime = cat(2, mTime, datenum(subStrs{2}, 'yyyy-mm-dd'));

    thisLine = fgetl(fid);
    subStrs = strsplit(thisLine, ' ');
    nProfiles = cat(2, nProfiles, str2double(subStrs{2}));

    fgetl(fid);

    fgetl(fid);
    fgetl(fid);
    thisLine = fgetl(fid);
    subStrs = strsplit(thisLine, ' ');
    meanErrFernald0297 = cat(2, meanErrFernald0297, str2double(subStrs{2}(1:(end - 2))));
    stdErrFernald0297 = cat(2, stdErrFernald0297, str2double(subStrs{4}(1:(end - 2))));
    outliersFernald0297 = cat(2, outliersFernald0297, str2double(subStrs{6}));

    thisLine = fgetl(fid);
    subStrs = strsplit(thisLine, ' ');
    meanErrXian0297 = cat(2, meanErrXian0297, str2double(subStrs{2}(1:(end - 2))));
    stdErrXian0297 = cat(2, stdErrXian0297, str2double(subStrs{4}(1:(end - 2))));
    outliersXian0297 = cat(2, outliersXian0297, str2double(subStrs{6}));

    fgetl(fid);

    thisLine = fgetl(fid);
    subStrs = strsplit(thisLine, ' ');
    meanErrLK0297 = cat(2, meanErrLK0297, str2double(subStrs{2}(1:(end - 2))));
    stdErrLK0297 = cat(2, stdErrLK0297, str2double(subStrs{4}(1:(end - 2))));
    outliersLK0297 = cat(2, outliersLK0297, str2double(subStrs{6}));

    fgetl(fid);
    fgetl(fid);
    thisLine = fgetl(fid);
    subStrs = strsplit(thisLine, ' ');
    meanErrFernald0998 = cat(2, meanErrFernald0998, str2double(subStrs{2}(1:(end - 2))));
    stdErrFernald0998 = cat(2, stdErrFernald0998, str2double(subStrs{4}(1:(end - 2))));
    outliersFernald0998 = cat(2, outliersFernald0998, str2double(subStrs{6}));

    thisLine = fgetl(fid);
    subStrs = strsplit(thisLine, ' ');
    meanErrXian0998 = cat(2, meanErrXian0998, str2double(subStrs{2}(1:(end - 2))));
    stdErrXian0998 = cat(2, stdErrXian0998, str2double(subStrs{4}(1:(end - 2))));
    outliersXian0998 = cat(2, outliersXian0998, str2double(subStrs{6}));

    fgetl(fid);

    thisLine = fgetl(fid);
    subStrs = strsplit(thisLine, ' ');
    meanErrLK0998 = cat(2, meanErrLK0998, str2double(subStrs{2}(1:(end - 2))));
    stdErrLK0998 = cat(2, stdErrLK0998, str2double(subStrs{4}(1:(end - 2))));
    outliersLK0998 = cat(2, outliersLK0998, str2double(subStrs{6}));

    fgetl(fid);
end

fclose(fid);

%% Display
figure('Position', [0, 20, 600, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;

p1 = plot(mTime, meanErrFernald0297, '-k', 'DisplayName', 'Fernald算法');
p2 = plot(mTime, meanErrXian0297, '-b', 'DisplayName', '大舜算法');
p4 = plot(mTime, meanErrLK0297, '-r', 'DisplayName', '蓝科光电');

plot(mTime, zeros(size(mTime)), 'LineStyle', '-.', 'color', 'cyan');

hold off;

xlim([min(mTime) - 1, max(mTime) + 1]);
ylim([-10, 10]);

xlabel('日期');
ylabel('偏差 (千米)');
title('能见度反演方法平均偏差（I0297）');

set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

datetick(gca, 'x', 'mm-dd', 'keeplimits', 'keepticks');

legend([p1, p2, p4], 'Location', 'NorthEast');

figure('Position', [0, 20, 600, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;

p1 = plot(mTime, outliersFernald0297, '-k', 'DisplayName', 'Fernald算法');
p2 = plot(mTime, outliersXian0297, '-b', 'DisplayName', '大舜算法');
p4 = plot(mTime, outliersLK0297, '-r', 'DisplayName', '蓝科光电');

plot(mTime, zeros(size(mTime)), 'LineStyle', '-.', 'color', 'cyan');

hold off;

xlim([min(mTime) - 1, max(mTime) + 1]);
ylim([0, 100]);

xlabel('日期');
ylabel('个数');
title('能见度反演异常廓线数（I0297）');

set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

datetick(gca, 'x', 'mm-dd', 'keeplimits', 'keepticks');

legend([p1, p2, p4], 'Location', 'NorthEast');


figure('Position', [0, 20, 600, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;

p1 = plot(mTime, stdErrFernald0297, '-k', 'DisplayName', 'Fernald算法');
p2 = plot(mTime, stdErrXian0297, '-b', 'DisplayName', '大舜算法');
p4 = plot(mTime, stdErrLK0297, '-r', 'DisplayName', '蓝科光电');

plot(mTime, zeros(size(mTime)), 'LineStyle', '-.', 'color', 'cyan');

hold off;

xlim([min(mTime) - 1, max(mTime) + 1]);
ylim([0, 10]);

xlabel('日期');
ylabel('偏差 (千米)');
title('能见度反演方法偏差标准差（I0297）');

set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

datetick(gca, 'x', 'mm-dd', 'keeplimits', 'keepticks');

legend([p1, p2, p4], 'Location', 'NorthEast');
