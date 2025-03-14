% 分析中科光博能见度反演结果
% 作者：殷振平
% 日期：2025-03-14

clc;
close all;

%% Parameter Definition
statsFile = 'D:\CMA-vis-lidar-assessment\highway-obs\quicklooks\GBQ-Cmp\vis_statistics_new.txt';
l1Folder = 'D:\CMA-vis-lidar-assessment\highway-obs\L1\new';   % 雷达产品文件目录
saveFolder = 'D:\CMA-vis-lidar-assessment\highway-obs\quicklooks\GBQ-Cmp';   % 输出结果目录
GBQL2Data = 'D:\CMA-vis-lidar-assessment\highway-obs\GBQ-ret';
visSensorFile = 'vis-sensor-data.mat';   % 前向散射能见度仪数据文件
visDiffThresh = 2e4;   % [m]
dist2_998 = 4.5;   % I0998站点与激光雷达距离（千米）
dist2_297 = 3.2;   % I0297站点与激光雷达距离（千米）
flagDisplayMeanStats = true;
flagDisplayDailyStats = false;
visible = 'off';

%% 读取前向散射能见度仪数据
load(visSensorFile);

%% 展示每日统计结果
if flagDisplayDailyStats

    % 打开统计文件句柄
    sFid = fopen(statsFile, 'w');

    % 读取每天的能见度数据
    gbqFiles = listfile(GBQL2Data, '\w*_vis_lidar_l2.mat', 1);
    matFiles = listfile(l1Folder, '\w*_exp.mat', 1);

    mTime = [];
    height = [];
    visMat_Holger = [];
    visMat_Fernald = [];
    visMat_GBQ = [];
    visMat_LK = [];
    for iFile = 1:length(matFiles)

        fprintf('Processing %s\n', matFiles{iFile});

        data = load(matFiles{iFile});
        height = data.height;
        thisMTime = data.mTime;
        thisVisMat_Holger = data.visMat_Holger;
        thisVisMat_Fernald = data.visMat_Fernald;
        thisVisMat_LK = data.visMat_LK;

        dataGBQ = load(gbqFiles{iFile});
        thisMTimeGBQ = dataGBQ.mTime;
        thisRangeGBQ = dataGBQ.range * 1e-3;
        tmpVisMat_GBQ = dataGBQ.vis_gbq;

        % interpolation
        [TIME, RANGE] = meshgrid(thisMTimeGBQ, thisRangeGBQ);
        thisVisMat_GBQ = interp2(TIME, RANGE, tmpVisMat_GBQ, thisMTime, height);

        %% I0297
        refIdx297 = find(height >= dist2_297, 1);
        visFernald297 = thisVisMat_Fernald(refIdx297, :);
        visGBQ297 = thisVisMat_GBQ(refIdx297, :);
        visLK297 = thisVisMat_LK(refIdx297, :);
        vis297Interp = interp1(vis297.mTime, vis297.vis, thisMTime);

        diffFernald297 = visFernald297 - vis297Interp;
        stdFernald297 = nanstd(diffFernald297((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)));
        meanFernald297 = nanmean(diffFernald297((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)));
        relStdFernald297 = nanstd(diffFernald297((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)) ./ vis297Interp((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)));
        relMeanFernald297 = nanmean(diffFernald297((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)) ./ vis297Interp((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)));
        nFernald297Outliers = sum(abs(diffFernald297) > visDiffThresh);

        diffGBQ297 = visGBQ297 - vis297Interp;
        stdGBQ297 = nanstd(diffGBQ297((abs(diffGBQ297) <= visDiffThresh) & (visGBQ297 > 0)));
        meanGBQ297 = nanmean(diffGBQ297((abs(diffGBQ297) <= visDiffThresh) & (visGBQ297 > 0)));
        relStdGBQ297 = nanstd(diffGBQ297((abs(diffGBQ297) <= visDiffThresh) & (visGBQ297 > 0)) ./ vis297Interp((abs(diffGBQ297) <= visDiffThresh) & (visGBQ297 > 0)));
        relMeanGBQ297 = nanmean(diffGBQ297((abs(diffGBQ297) <= visDiffThresh) & (visGBQ297 > 0)) ./ vis297Interp((abs(diffGBQ297) <= visDiffThresh) & (visGBQ297 > 0)));
        nGBQ297Outliers = sum(abs(diffGBQ297) > visDiffThresh);

        diffLK297 = visLK297 - vis297Interp;
        stdLK297 = nanstd(diffLK297((abs(diffLK297) <= visDiffThresh) & (visLK297 > 0)));
        meanLK297 = nanmean(diffLK297((abs(diffLK297) <= visDiffThresh) & (visLK297 > 0)));
        relStdLK297 = nanstd(diffLK297((abs(diffLK297) <= visDiffThresh) & (visLK297 > 0)) ./ vis297Interp((abs(diffLK297) <= visDiffThresh) & (visLK297 > 0)));
        relMeanLK297 = nanmean(diffLK297((abs(diffLK297) <= visDiffThresh) & (visLK297 > 0)) ./ vis297Interp((abs(diffLK297) <= visDiffThresh) & (visLK297 > 0)));
        nLK297Outliers = sum(abs(diffLK297) > visDiffThresh);

        % write statistics
        fprintf(sFid, 'Date: %s\n', datestr(thisMTime(1), 'yyyy-mm-dd'));
        fprintf(sFid, 'Total Profiles: %d\n\n', length(thisMTime));
        fprintf(sFid, 'I0297\n');
        fprintf(sFid, '方法: 平均偏差 平均相对偏差 标准差 平均相对偏差标准差 异常点个数(偏差超过%6.2fkm)\n', visDiffThresh * 1e-3);
        fprintf(sFid, 'Fernald: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanFernald297 * 1e-3, relMeanFernald297 * 100, stdFernald297 * 1e-3, relStdFernald297 * 100, nFernald297Outliers);
        fprintf(sFid, 'GBQ: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanGBQ297 * 1e-3, relMeanGBQ297 * 100, stdGBQ297 * 1e-3, relStdGBQ297 * 100, nGBQ297Outliers);
        fprintf(sFid, 'LK: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanLK297 * 1e-3, relMeanLK297 * 100, stdLK297 * 1e-3, relStdLK297 * 100, nLK297Outliers);

        % Display Statiscs
        figure('Position', [0, 30, 550, 500], 'Units', 'Pixels', 'Color', 'w', 'visible', visible);

        subplot('Position', [0.15, 0.6, 0.8, 0.35], 'Units', 'Normalized');

        hold on;
        p1 = plot(thisMTime, thisVisMat_Fernald(refIdx297, :) * 1e-3, '-k', 'DisplayName', 'Fernald算法');
        p2 = plot(thisMTime, thisVisMat_GBQ(refIdx297, :) * 1e-3, '-b', 'DisplayName', '光博算法');
        p4 = plot(thisMTime, thisVisMat_LK(refIdx297, :) * 1e-3, '-r', 'DisplayName', '蓝科光电');
        p5 = scatter(vis297.mTime, vis297.vis * 1e-3, 10, 'Marker', 's', 'MarkerEdgeColor', 'cyan', 'MarkerFaceColor', 'cyan', 'DisplayName', 'I0297');
        hold off;

        xlim([min(thisMTime), max(thisMTime)]);
        ylim([0, 50]);

        % xlabel('时间');
        ylabel('能见度 (千米)');
        title(sprintf('%s 能见度对比 (I0297)', datestr(thisMTime(1), 'yyyy-mm-dd')));

        set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

        datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

        legend([p1, p2, p4, p5], 'Location', 'NorthEast');

        subplot('Position', [0.15, 0.1, 0.8, 0.25], 'Units', 'Normalized');

        hold on;

        plot(thisMTime, diffFernald297 * 1e-3, '-k');
        plot(thisMTime, diffGBQ297 * 1e-3, '-b');
        plot(thisMTime, diffLK297 * 1e-3, '-r');

        plot(thisMTime, zeros(size(thisMTime)), 'LineStyle', '-.', 'color', 'cyan');

        hold off;

        xlim([min(thisMTime), max(thisMTime)]);
        ylim([-10, 10]);

        xlabel('时间');
        ylabel('能见度偏差 (千米)');

        set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

        datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

        text(0, 1.24, sprintf('(Fernald 算法) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(光博算法) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(蓝科光电) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n', ...
            meanFernald297 * 1e-3, relMeanFernald297 * 100, stdFernald297 * 1e-3, relStdFernald297 * 100, ...
            meanGBQ297 * 1e-3, relMeanGBQ297 * 100, stdGBQ297 * 1e-3, relStdGBQ297 * 100, ...
            meanLK297 * 1e-3, relMeanLK297 * 100, stdLK297 * 1e-3, relStdLK297 * 100), 'Units', 'Normalized', 'FontSize', 10);

        export_fig(gcf, fullfile(saveFolder, sprintf('%s_vis_comparison_297.png', datestr(thisMTime(1), 'yyyy-mm-dd'))), '-r300');

        %% Concatnate data
        mTime = cat(2, mTime, thisMTime);
        visMat_GBQ = cat(2, visMat_GBQ, thisVisMat_GBQ);
        visMat_Fernald = cat(2, visMat_Fernald, thisVisMat_Fernald);
        visMat_LK = cat(2, visMat_LK, thisVisMat_LK);

    end

    fclose(sFid);
end

%% 展示天平均统计结果
if flagDisplayMeanStats
    mTime = [];
    nProfiles = [];
    meanErrFernald0297 = [];
    stdErrFernald0297 = [];
    outliersFernald0297 = [];
    meanErrLK0297 = [];
    stdErrLK0297 = [];
    outliersLK0297 = [];
    meanErrGBQ0297 = [];
    stdErrGBQ0297 = [];
    outliersGBQ0297 = [];

    %% Read data
    fid = fopen(statsFile, 'r');

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
        meanErrGBQ0297 = cat(2, meanErrGBQ0297, str2double(subStrs{2}(1:(end - 2))));
        stdErrGBQ0297 = cat(2, stdErrGBQ0297, str2double(subStrs{4}(1:(end - 2))));
        outliersGBQ0297 = cat(2, outliersGBQ0297, str2double(subStrs{6}));

        thisLine = fgetl(fid);
        subStrs = strsplit(thisLine, ' ');
        meanErrLK0297 = cat(2, meanErrLK0297, str2double(subStrs{2}(1:(end - 2))));
        stdErrLK0297 = cat(2, stdErrLK0297, str2double(subStrs{4}(1:(end - 2))));
        outliersLK0297 = cat(2, outliersLK0297, str2double(subStrs{6}));
    end

    fclose(fid);

    %% Display
    figure('Position', [0, 20, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'visible', visible);

    hold on;

    p1 = plot(mTime, meanErrFernald0297, '-k', 'Marker', '*', 'DisplayName', 'Fernald算法');
    p2 = plot(mTime, meanErrGBQ0297, '-b', 'Marker', '*', 'DisplayName', '光博算法');
    p4 = plot(mTime, meanErrLK0297, '-r', 'Marker', '*', 'DisplayName', '蓝科光电');

    plot(mTime, zeros(size(mTime)), 'LineStyle', '-.', 'color', 'cyan');

    hold off;

    xlim([min(mTime) - 1, max(mTime) + 1]);
    ylim([-10, 10]);

    xlabel('日期');
    ylabel('偏差 (千米)');
    title('能见度反演方法平均偏差（I0297）');
    text(0.02, 0.89, sprintf('平均误差\nFernald: %5.2fkm; 光博:%5.2fkm; 蓝科:%5.2fkm', nanmean(meanErrFernald0297), nanmean(meanErrGBQ0297), nanmean(meanErrLK0297)), 'FontSize', 11, 'Units', 'Normalized');

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

    datetick(gca, 'x', 'mm-dd', 'keeplimits', 'keepticks');

    legend([p1, p2, p4], 'Location', 'NorthEast');

    export_fig(gcf, fullfile(saveFolder, sprintf('Overview_vis_comparison_mean_err_297.png')), '-r300');

    figure('Position', [0, 20, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'visible', visible);

    hold on;

    p1 = plot(mTime, outliersFernald0297, '-k', 'Marker', '*', 'DisplayName', 'Fernald算法');
    p2 = plot(mTime, outliersGBQ0297, '-b', 'Marker', '*', 'DisplayName', '光博算法');
    p4 = plot(mTime, outliersLK0297, '-r', 'Marker', '*', 'DisplayName', '蓝科光电');

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

    export_fig(gcf, fullfile(saveFolder, sprintf('Overview_vis_comparison_outliers_297.png')), '-r300');

    figure('Position', [0, 20, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'visible', visible);

    hold on;

    p1 = plot(mTime, stdErrFernald0297, '-k', 'Marker', '*', 'DisplayName', 'Fernald算法');
    p2 = plot(mTime, stdErrGBQ0297, '-b', 'Marker', '*', 'DisplayName', '光博算法');
    p4 = plot(mTime, stdErrLK0297, '-r', 'Marker', '*', 'DisplayName', '蓝科光电');

    % plot(mTime, zeros(size(mTime)), 'LineStyle', '-.', 'color', 'cyan');

    hold off;

    xlim([min(mTime) - 1, max(mTime) + 1]);
    ylim([0, 10]);

    xlabel('日期');
    ylabel('偏差 (千米)');
    title('能见度反演方法偏差标准差（I0297）');

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);
    text(0.02, 0.86, sprintf('标准偏差\nFernald: %5.2fkm; 光博:%5.2fkm; 蓝科:%5.2fkm', nanmean(stdErrFernald0297), nanmean(stdErrGBQ0297), nanmean(stdErrLK0297)), 'FontSize', 11, 'Units', 'Normalized');

    datetick(gca, 'x', 'mm-dd', 'keeplimits', 'keepticks');

    legend([p1, p2, p4], 'Location', 'NorthEast');

    export_fig(gcf, fullfile(saveFolder, sprintf('Overview_vis_comparison_std_err_297.png')), '-r300');

end