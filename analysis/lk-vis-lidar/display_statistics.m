% 分析能见度反演结果
% 日期：2025-02-05

clc;
close all;

%% Parameter Definition
statsFile = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\quicklooks\vis_statistics_new.txt';
l1Folder = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\L1\new';   % 雷达产品文件目录
saveFolder = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\quicklooks';   % 输出结果目录
visSensorFile = 'vis-sensor-data.mat';   % 前向散射能见度仪数据文件
visDiffThresh = 2e4;   % [m]
dist2_998 = 4.5;   % I0998站点与激光雷达距离（千米）
dist2_297 = 3.2;   % I0297站点与激光雷达距离（千米）
flagDisplayMeanStats = true;
flagDisplayDailyStats = true;
visible = 'off';

%% 读取前向散射能见度仪数据
load(visSensorFile);

%% 展示每日统计结果
if flagDisplayDailyStats

    % 打开统计文件句柄
    sFid = fopen(statsFile, 'w');

    % 读取每天的能见度数据
    matFiles = listfile(l1Folder, '\w*_exp.mat', 1);

    mTime = [];
    height = [];
    visMat_Xian = [];
    visMat_Holger = [];
    visMat_Fernald = [];
    visMat_LK = [];
    for iFile = 1:length(matFiles)

        fprintf('Processing %s\n', matFiles{iFile});

        data = load(matFiles{iFile});
        height = data.height;
        thisMTime = data.mTime;
        thisVisMat_Xian = data.visMat_Xian;
        thisVisMat_Holger = data.visMat_Holger;
        thisVisMat_Fernald = data.visMat_Fernald;
        thisVisMat_LK = data.visMat_LK;

        %% I0297
        refIdx297 = find(height >= dist2_297, 1);
        visFernald297 = thisVisMat_Fernald(refIdx297, :);
        visXian297 = thisVisMat_Xian(refIdx297, :);
        visHolger297 = thisVisMat_Holger(refIdx297, :);
        visLK297 = thisVisMat_LK(refIdx297, :);
        vis297Interp = interp1(vis297.mTime, vis297.vis, thisMTime);

        diffFernald297 = visFernald297 - vis297Interp;
        stdFernald297 = nanstd(diffFernald297((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)));
        meanFernald297 = nanmean(diffFernald297((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)));
        relStdFernald297 = nanstd(diffFernald297((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)) ./ vis297Interp((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)));
        relMeanFernald297 = nanmean(diffFernald297((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)) ./ vis297Interp((abs(diffFernald297) <= visDiffThresh) & (visFernald297 > 0)));
        nFernald297Outliers = sum(abs(diffFernald297) > visDiffThresh);

        diffXian297 = visXian297 - vis297Interp;
        stdXian297 = nanstd(diffXian297((abs(diffXian297) <= visDiffThresh) & (visXian297 > 0)));
        meanXian297 = nanmean(diffXian297((abs(diffXian297) <= visDiffThresh) & (visXian297 > 0)));
        relStdXian297 = nanstd(diffXian297((abs(diffXian297) <= visDiffThresh) & (visXian297 > 0)) ./ vis297Interp((abs(diffXian297) <= visDiffThresh) & (visXian297 > 0)));
        relMeanXian297 = nanmean(diffXian297((abs(diffXian297) <= visDiffThresh) & (visXian297 > 0)) ./ vis297Interp((abs(diffXian297) <= visDiffThresh) & (visXian297 > 0)));
        nXian297Outliers = sum(abs(diffXian297) > visDiffThresh);

        diffHolger297 = visHolger297 - vis297Interp;
        stdHolger297 = nanstd(diffHolger297((abs(diffHolger297) <= visDiffThresh) & (visHolger297 > 0)));
        meanHolger297 = nanmean(diffHolger297((abs(diffHolger297) <= visDiffThresh) & (visHolger297 > 0)));
        relStdHolger297 = nanstd(diffHolger297((abs(diffHolger297) <= visDiffThresh) & (visHolger297 > 0)) ./ vis297Interp((abs(diffHolger297) <= visDiffThresh) & (visHolger297 > 0)));
        relMeanHolger297 = nanmean(diffHolger297((abs(diffHolger297) <= visDiffThresh) & (visHolger297 > 0)) ./ vis297Interp((abs(diffHolger297) <= visDiffThresh) & (visHolger297 > 0)));
        nHolger297Outliers = sum(abs(diffHolger297) > visDiffThresh);

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
        fprintf(sFid, 'Xian: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanXian297 * 1e-3, relMeanXian297 * 100, stdXian297 * 1e-3, relStdXian297 * 100, nXian297Outliers);
        fprintf(sFid, 'Holger: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanHolger297 * 1e-3, relMeanHolger297 * 100, stdHolger297 * 1e-3, relStdHolger297 * 100, nHolger297Outliers);
        fprintf(sFid, 'LK: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanLK297 * 1e-3, relMeanLK297 * 100, stdLK297 * 1e-3, relStdLK297 * 100, nLK297Outliers);

        % Display Statiscs
        figure('Position', [0, 30, 550, 500], 'Units', 'Pixels', 'Color', 'w', 'visible', visible);

        subplot('Position', [0.15, 0.6, 0.8, 0.35], 'Units', 'Normalized');

        hold on;
        p1 = plot(thisMTime, thisVisMat_Fernald(refIdx297, :) * 1e-3, '-k', 'DisplayName', 'Fernald算法');
        p2 = plot(thisMTime, thisVisMat_Xian(refIdx297, :) * 1e-3, '-b', 'DisplayName', '大舜算法');
        p3 = plot(thisMTime, thisVisMat_Holger(refIdx297, :) * 1e-3, '-g', 'DisplayName', '前向迭代');
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

        legend([p1, p2, p3, p4, p5], 'Location', 'NorthEast');

        subplot('Position', [0.15, 0.1, 0.8, 0.25], 'Units', 'Normalized');

        hold on;

        plot(thisMTime, diffFernald297 * 1e-3, '-k');
        plot(thisMTime, diffXian297 * 1e-3, '-b');
        plot(thisMTime, diffHolger297 * 1e-3, '-g');
        plot(thisMTime, diffLK297 * 1e-3, '-r');

        plot(thisMTime, zeros(size(thisMTime)), 'LineStyle', '-.', 'color', 'cyan');

        hold off;

        xlim([min(thisMTime), max(thisMTime)]);
        ylim([-10, 10]);

        xlabel('时间');
        ylabel('能见度偏差 (千米)');

        set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

        datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

        text(0, 1.24, sprintf('(Fernald 算法) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(大舜算法) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(前向迭代) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(蓝科光电) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n', ...
            meanFernald297 * 1e-3, relMeanFernald297 * 100, stdFernald297 * 1e-3, relStdFernald297 * 100, ...
            meanXian297 * 1e-3, relMeanXian297 * 100, stdXian297 * 1e-3, relStdXian297 * 100, ...
            meanHolger297 * 1e-3, relMeanHolger297 * 100, stdHolger297 * 1e-3, relStdHolger297 * 100, ...
            meanLK297 * 1e-3, relMeanLK297 * 100, stdLK297 * 1e-3, relStdLK297 * 100), 'Units', 'Normalized', 'FontSize', 10);

        export_fig(gcf, fullfile(saveFolder, sprintf('%s_vis_comparison_297.png', datestr(thisMTime(1), 'yyyy-mm-dd'))), '-r300');

        % % laser temperature
        % figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w', 'visible', visible);

        % hold on;
        % plot(thisMTime, laserTemp, '-k');
        % hold off;

        % xlim([min(thisMTime), max(thisMTime)]);
        % ylim([0, 50]);

        % xlabel('时间');
        % ylabel('温度（度）');
        % title(sprintf('%s 激光头温度', datestr(thisMTime(1), 'yyyy-mm-dd')));

        % set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

        % datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

        % export_fig(gcf, fullfile(saveFolder, sprintf('%s_laser_temperature.png', datestr(thisMTime(1), 'yyyy-mm-dd'))), '-r300');

        %% I0998
        refIdx998 = find(height >= dist2_998, 1);
        visFernald998 = thisVisMat_Fernald(refIdx998, :);
        visXian998 = thisVisMat_Xian(refIdx998, :);
        visHolger998 = thisVisMat_Holger(refIdx998, :);
        visLK998 = thisVisMat_LK(refIdx998, :);
        vis998Interp = interp1(vis998.mTime, vis998.vis, thisMTime);

        diffFernald998 = visFernald998 - vis998Interp;
        stdFernald998 = nanstd(diffFernald998((abs(diffFernald998) <= visDiffThresh) & (visFernald998 > 0)));
        meanFernald998 = nanmean(diffFernald998((abs(diffFernald998) <= visDiffThresh) & (visFernald998 > 0)));
        relStdFernald998 = nanstd(diffFernald998((abs(diffFernald998) <= visDiffThresh) & (visFernald998 > 0)) ./ vis998Interp((abs(diffFernald998) <= visDiffThresh) & (visFernald998 > 0)));
        relMeanFernald998 = nanmean(diffFernald998((abs(diffFernald998) <= visDiffThresh) & (visFernald998 > 0)) ./ vis998Interp((abs(diffFernald998) <= visDiffThresh) & (visFernald998 > 0)));
        nFernald998Outliers = sum(abs(diffFernald998) > visDiffThresh);

        diffXian998 = visXian998 - vis998Interp;
        stdXian998 = nanstd(diffXian998((abs(diffXian998) <= visDiffThresh) & (visXian998 > 0)));
        meanXian998 = nanmean(diffXian998((abs(diffXian998) <= visDiffThresh) & (visXian998 > 0)));
        relStdXian998 = nanstd(diffXian998((abs(diffXian998) <= visDiffThresh) & (visXian998 > 0)) ./ vis998Interp((abs(diffXian998) <= visDiffThresh) & (visXian998 > 0)));
        relMeanXian998 = nanmean(diffXian998((abs(diffXian998) <= visDiffThresh) & (visXian998 > 0)) ./ vis998Interp((abs(diffXian998) <= visDiffThresh) & (visXian998 > 0)));
        nXian998Outliers = sum(abs(diffXian998) > visDiffThresh);

        diffHolger998 = visHolger998 - vis998Interp;
        stdHolger998 = nanstd(diffHolger998((abs(diffHolger998) <= visDiffThresh) & (visHolger998 > 0)));
        meanHolger998 = nanmean(diffHolger998((abs(diffHolger998) <= visDiffThresh) & (visHolger998 > 0)));
        relStdHolger998 = nanstd(diffHolger998((abs(diffHolger998) <= visDiffThresh) & (visHolger998 > 0)) ./ vis998Interp((abs(diffHolger998) <= visDiffThresh) & (visHolger998 > 0)));
        relMeanHolger998 = nanmean(diffHolger998((abs(diffHolger998) <= visDiffThresh) & (visHolger998 > 0)) ./ vis998Interp((abs(diffHolger998) <= visDiffThresh) & (visHolger998 > 0)));
        nHolger998Outliers = sum(abs(diffHolger998) > visDiffThresh);

        diffLK998 = visLK998 - vis998Interp;
        stdLK998 = nanstd(diffLK998((abs(diffLK998) <= visDiffThresh) & (visLK998 > 0)));
        meanLK998 = nanmean(diffLK998((abs(diffLK998) <= visDiffThresh) & (visLK998 > 0)));
        relStdLK998 = nanstd(diffLK998((abs(diffLK998) <= visDiffThresh) & (visLK998 > 0)) ./ vis998Interp((abs(diffLK998) <= visDiffThresh) & (visLK998 > 0)));
        relMeanLK998 = nanmean(diffLK998((abs(diffLK998) <= visDiffThresh) & (visLK998 > 0)) ./ vis998Interp((abs(diffLK998) <= visDiffThresh) & (visLK998 > 0)));
        nLK998Outliers = sum(abs(diffLK998) > visDiffThresh);

        % write statistics
        fprintf(sFid, 'I0998\n');
        fprintf(sFid, '方法: 平均偏差 平均相对偏差 标准差 平均相对偏差标准差 异常点个数(偏差超过%6.2fkm)\n', visDiffThresh * 1e-3);
        fprintf(sFid, 'Fernald: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanFernald998 * 1e-3, relMeanFernald998 * 100 , stdFernald998 * 1e-3, relStdFernald998 * 100, nFernald998Outliers);
        fprintf(sFid, 'Xian: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanXian998 * 1e-3, relMeanXian998 * 100, stdXian998 * 1e-3, relStdXian998 * 100, nXian998Outliers);
        fprintf(sFid, 'Holger: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanHolger998 * 1e-3, relMeanHolger998 * 100, stdHolger998 * 1e-3, relStdHolger998 * 100, nHolger998Outliers);
        fprintf(sFid, 'LK: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n\n', meanLK998 * 1e-3, relMeanLK998 * 100, stdLK998 * 1e-3, relStdLK998 * 100, nLK998Outliers);

        % Display Statiscs
        figure('Position', [0, 30, 550, 500], 'Units', 'Pixels', 'Color', 'w', 'visible', visible);

        subplot('Position', [0.15, 0.6, 0.8, 0.35], 'Units', 'Normalized');

        hold on;
        p1 = plot(thisMTime, thisVisMat_Fernald(refIdx998, :) * 1e-3, '-k', 'DisplayName', 'Fernald算法');
        p2 = plot(thisMTime, thisVisMat_Xian(refIdx998, :) * 1e-3, '-b', 'DisplayName', '大舜算法');
        p3 = plot(thisMTime, thisVisMat_Holger(refIdx998, :) * 1e-3, '-g', 'DisplayName', '前向迭代');
        p4 = plot(thisMTime, thisVisMat_LK(refIdx998, :) * 1e-3, '-r', 'DisplayName', '蓝科光电');
        p5 = scatter(vis998.mTime, vis998.vis * 1e-3, 10, 'Marker', 's', 'MarkerEdgeColor', 'cyan', 'MarkerFaceColor', 'cyan', 'DisplayName', 'I0998');
        hold off;

        xlim([min(thisMTime), max(thisMTime)]);
        ylim([0, 50]);

        % xlabel('时间');
        ylabel('能见度 (千米)');
        title(sprintf('%s 能见度对比 (I0998)', datestr(thisMTime(1), 'yyyy-mm-dd')));

        set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

        datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

        legend([p1, p2, p3, p4, p5], 'Location', 'NorthEast');

        subplot('Position', [0.15, 0.1, 0.8, 0.25], 'Units', 'Normalized');

        hold on;

        plot(thisMTime, diffFernald998 * 1e-3, '-k');
        plot(thisMTime, diffXian998 * 1e-3, '-b');
        plot(thisMTime, diffHolger998 * 1e-3, '-g');
        plot(thisMTime, diffLK998 * 1e-3, '-r');

        plot(thisMTime, zeros(size(thisMTime)), 'LineStyle', '-.', 'color', 'cyan');

        hold off;

        xlim([min(thisMTime), max(thisMTime)]);
        ylim([-10, 10]);

        xlabel('时间');
        ylabel('能见度偏差 (千米)');

        set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

        datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

        text(0, 1.24, sprintf('(Fernald 算法) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(大舜算法) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(前向迭代) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(蓝科光电) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n', ...
            meanFernald998 * 1e-3, relMeanFernald998 * 100, stdFernald998 * 1e-3, relStdFernald998 * 100, ...
            meanXian998 * 1e-3, relMeanXian998 * 100, stdXian998 * 1e-3, relStdXian998 * 100, ...
            meanHolger998 * 1e-3, relMeanHolger998 * 100, stdHolger998 * 1e-3, relStdHolger998 * 100, ...
            meanLK998 * 1e-3, relMeanLK998 * 100, stdLK998 * 1e-3, relStdLK998 * 100), 'Units', 'Normalized', 'FontSize', 10);

        export_fig(gcf, fullfile(saveFolder, sprintf('%s_vis_comparison_998.png', datestr(thisMTime(1), 'yyyy-mm-dd'))), '-r300');

        %% Concatnate data
        mTime = cat(2, mTime, thisMTime);
        visMat_Xian = cat(2, visMat_Xian, thisVisMat_Xian);
        visMat_Holger = cat(2, visMat_Holger, thisVisMat_Holger);
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
    figure('Position', [0, 20, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'visible', visible);

    hold on;

    p1 = plot(mTime, meanErrFernald0297, '-k', 'Marker', '*', 'DisplayName', 'Fernald算法');
    p2 = plot(mTime, meanErrXian0297, '-b', 'Marker', '*', 'DisplayName', '大舜算法');
    p4 = plot(mTime, meanErrLK0297, '-r', 'Marker', '*', 'DisplayName', '蓝科光电');

    plot(mTime, zeros(size(mTime)), 'LineStyle', '-.', 'color', 'cyan');

    hold off;

    xlim([min(mTime) - 1, max(mTime) + 1]);
    ylim([-10, 10]);

    xlabel('日期');
    ylabel('偏差 (千米)');
    title('能见度反演方法平均偏差（I0297）');
    text(0.02, 0.89, sprintf('平均误差\nFernald: %5.2fkm; 大舜:%5.2fkm; 蓝科:%5.2fkm', nanmean(meanErrFernald0297), nanmean(meanErrXian0297), nanmean(meanErrLK0297)), 'FontSize', 11, 'Units', 'Normalized');

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

    datetick(gca, 'x', 'mm-dd', 'keeplimits', 'keepticks');

    legend([p1, p2, p4], 'Location', 'NorthEast');

    export_fig(gcf, fullfile(saveFolder, sprintf('Overview_vis_comparison_mean_err_297.png')), '-r300');

    figure('Position', [0, 20, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'visible', visible);

    hold on;

    p1 = plot(mTime, outliersFernald0297, '-k', 'Marker', '*', 'DisplayName', 'Fernald算法');
    p2 = plot(mTime, outliersXian0297, '-b', 'Marker', '*', 'DisplayName', '大舜算法');
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
    p2 = plot(mTime, stdErrXian0297, '-b', 'Marker', '*', 'DisplayName', '大舜算法');
    p4 = plot(mTime, stdErrLK0297, '-r', 'Marker', '*', 'DisplayName', '蓝科光电');

    % plot(mTime, zeros(size(mTime)), 'LineStyle', '-.', 'color', 'cyan');

    hold off;

    xlim([min(mTime) - 1, max(mTime) + 1]);
    ylim([0, 10]);

    xlabel('日期');
    ylabel('偏差 (千米)');
    title('能见度反演方法偏差标准差（I0297）');

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);
    text(0.02, 0.86, sprintf('标准偏差\nFernald: %5.2fkm; 大舜:%5.2fkm; 蓝科:%5.2fkm', nanmean(stdErrFernald0297), nanmean(stdErrXian0297), nanmean(stdErrLK0297)), 'FontSize', 11, 'Units', 'Normalized');

    datetick(gca, 'x', 'mm-dd', 'keeplimits', 'keepticks');

    legend([p1, p2, p4], 'Location', 'NorthEast');

    export_fig(gcf, fullfile(saveFolder, sprintf('Overview_vis_comparison_std_err_297.png')), '-r300');

end