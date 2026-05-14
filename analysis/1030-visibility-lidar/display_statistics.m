% 分析能见度反演结果
% 日期：2025-02-05

clc;
close all;

%% Parameter Definition
statsFile = 'G:\backup\vis-lidar\20260428-tianjin-evaluation\前散仪对比数据\VIS1\vis_statistics_new.txt';
l1Folder = 'G:\backup\vis-lidar\20260428-tianjin-evaluation\前散仪对比数据\VIS1\L1';   % 雷达产品文件目录
saveFolder = 'H:\research\vislidar-intercomparison\tianjing\VIS1';   % 输出结果目录
visSensorFile = 'vis-sensor-data.mat';   % 前向散射能见度仪数据文件
visDiffThresh = 2e4;   % [m]
dist_vis_lidar = 0.5;   % 前向散射能见度仪与激光雷达距离（千米）
flagDisplayMeanStats = true;
flagDisplayDailyStats = true;
visible = 'on';

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
    for iFile = 1:length(matFiles)

        fprintf('Processing %s\n', matFiles{iFile});

        data = load(matFiles{iFile});
        height = data.range;
        thisMTime = data.mTime;
        thisVisMat_Xian = data.visMat_Xian;
        thisVisMat_Holger = data.visMat_Holger;
        thisVisMat_Fernald = data.visMat_Fernald;

        refIdx = find(height * 1e-3 >= dist_vis_lidar, 1);
        visFernald = thisVisMat_Fernald(refIdx, :);
        visXian = thisVisMat_Xian(refIdx, :);
        visHolger = thisVisMat_Holger(refIdx, :);
        visDataInterp = interp1(visData.mTime, visData.vis, thisMTime);

        diffFernald = visFernald - visDataInterp;
        stdFernald = nanstd(diffFernald((abs(diffFernald) <= visDiffThresh) & (visFernald > 0)));
        meanFernald = nanmean(diffFernald((abs(diffFernald) <= visDiffThresh) & (visFernald > 0)));
        relStdFernald = nanstd(diffFernald((abs(diffFernald) <= visDiffThresh) & (visFernald > 0)) ./ visDataInterp((abs(diffFernald) <= visDiffThresh) & (visFernald > 0)));
        relMeanFernald = nanmean(diffFernald((abs(diffFernald) <= visDiffThresh) & (visFernald > 0)) ./ visDataInterp((abs(diffFernald) <= visDiffThresh) & (visFernald > 0)));
        nFernaldOutliers = sum(abs(diffFernald) > visDiffThresh);

        diffXian = visXian - visDataInterp;
        stdXian = nanstd(diffXian((abs(diffXian) <= visDiffThresh) & (visXian > 0)));
        meanXian = nanmean(diffXian((abs(diffXian) <= visDiffThresh) & (visXian > 0)));
        relStdXian = nanstd(diffXian((abs(diffXian) <= visDiffThresh) & (visXian > 0)) ./ visDataInterp((abs(diffXian) <= visDiffThresh) & (visXian > 0)));
        relMeanXian = nanmean(diffXian((abs(diffXian) <= visDiffThresh) & (visXian > 0)) ./ visDataInterp((abs(diffXian) <= visDiffThresh) & (visXian > 0)));
        nXianOutliers = sum(abs(diffXian) > visDiffThresh);

        diffHolger = visHolger - visDataInterp;
        stdHolger = nanstd(diffHolger((abs(diffHolger) <= visDiffThresh) & (visHolger > 0)));
        meanHolger = nanmean(diffHolger((abs(diffHolger) <= visDiffThresh) & (visHolger > 0)));
        relStdHolger = nanstd(diffHolger((abs(diffHolger) <= visDiffThresh) & (visHolger > 0)) ./ visDataInterp((abs(diffHolger) <= visDiffThresh) & (visHolger > 0)));
        relMeanHolger = nanmean(diffHolger((abs(diffHolger) <= visDiffThresh) & (visHolger > 0)) ./ visDataInterp((abs(diffHolger) <= visDiffThresh) & (visHolger > 0)));
        nHolgerOutliers = sum(abs(diffHolger) > visDiffThresh);

        % write statistics
        fprintf(sFid, 'Date: %s\n', datestr(thisMTime(1), 'yyyy-mm-dd'));
        fprintf(sFid, 'Total Profiles: %d\n\n', length(thisMTime));
        fprintf(sFid, 'A6001\n');
        fprintf(sFid, '方法: 平均偏差 平均相对偏差 标准差 平均相对偏差标准差 异常点个数(偏差超过%6.2fkm)\n', visDiffThresh * 1e-3);
        fprintf(sFid, 'Fernald: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanFernald * 1e-3, relMeanFernald * 100, stdFernald * 1e-3, relStdFernald * 100, nFernaldOutliers);
        fprintf(sFid, 'Xian: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanXian * 1e-3, relMeanXian * 100, stdXian * 1e-3, relStdXian * 100, nXianOutliers);
        fprintf(sFid, 'Holger: %5.2fkm %6.2f%% %5.2fkm %6.2f%% %d\n', meanHolger * 1e-3, relMeanHolger * 100, stdHolger * 1e-3, relStdHolger * 100, nHolgerOutliers);

        % Display Statiscs
        figure('Position', [0, 30, 550, 500], 'Units', 'Pixels', 'Color', 'w', 'visible', visible);

        subplot('Position', [0.15, 0.6, 0.8, 0.35], 'Units', 'Normalized');

        hold on;
        p1 = plot(thisMTime, thisVisMat_Fernald(refIdx, :) * 1e-3, '-k', 'DisplayName', 'Fernald算法');
        p2 = plot(thisMTime, thisVisMat_Xian(refIdx, :) * 1e-3, '-b', 'DisplayName', '大舜算法');
        p3 = plot(thisMTime, thisVisMat_Holger(refIdx, :) * 1e-3, '-g', 'DisplayName', '前向迭代');
        p5 = scatter(visData.mTime, visData.vis * 1e-3, 10, 'Marker', 's', 'MarkerEdgeColor', 'cyan', 'MarkerFaceColor', 'cyan', 'DisplayName', '前散.');
        hold off;

        xlim([min(thisMTime), max(thisMTime)]);
        ylim([0, 50]);

        % xlabel('时间');
        ylabel('能见度 (千米)');
        title(sprintf('%s 能见度对比 (I0)', datestr(thisMTime(1), 'yyyy-mm-dd')));

        set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

        datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

        legend([p1, p2, p3, p5], 'Location', 'NorthEast');

        subplot('Position', [0.15, 0.1, 0.8, 0.25], 'Units', 'Normalized');

        hold on;

        plot(thisMTime, diffFernald * 1e-3, '-k');
        plot(thisMTime, diffXian * 1e-3, '-b');
        plot(thisMTime, diffHolger * 1e-3, '-g');

        plot(thisMTime, zeros(size(thisMTime)), 'LineStyle', '-.', 'color', 'cyan');

        hold off;

        xlim([min(thisMTime), max(thisMTime)]);
        ylim([-10, 10]);

        xlabel('时间');
        ylabel('能见度偏差 (千米)');

        set(gca, 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);

        datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

        text(0, 1.24, sprintf('(Fernald 算法) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(大舜算法) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n(前向迭代) 平均偏差: %5.2fkm (%6.2f%%); 标准差: %5.2fkm (%6.2f%%);\n', ...
            meanFernald * 1e-3, relMeanFernald * 100, stdFernald * 1e-3, relStdFernald * 100, ...
            meanXian * 1e-3, relMeanXian * 100, stdXian * 1e-3, relStdXian * 100, ...
            meanHolger * 1e-3, relMeanHolger * 100, stdHolger * 1e-3, relStdHolger * 100), ...
            'Units', 'Normalized', 'FontSize', 10);

        export_fig(gcf, fullfile(saveFolder, sprintf('%s_vis_comparison.png', datestr(thisMTime(1), 'yyyy-mm-dd'))), '-r300');

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

    end
    fclose(sFid);
end