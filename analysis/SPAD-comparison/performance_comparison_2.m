% 对比评估麓邦光电和Excelitas公司单光子探测器的性能
% 2025-09-16

%% Parameter Definition
dataFolder = 'C:\Users\zhenp\Documents\Industry\麓邦光电\LBTEK_APD';
outputFolder = 'C:\Users\zhenp\Documents\Industry\麓邦光电\LBTEK_APD\results';

%% Dark Count Comparison
subDataPath = fullfile(dataFolder, '1_测暗计数_Threshold_one_point_two_five');
files = listfile(subDataPath, '.*dat', 1);

for iFile = 1:length(files)
    lidarData = readALADat(files{iFile}, 'nMaxBin', 1900);
    lbtekData = lidarData.rawSignal(:, 1);
    exceData = lidarData.rawSignal(:, 2);

    sumDarkCountLBTEK = sum(lbtekData);
    sumDarkCountExcelitas = sum(exceData);

    figure('Position', [0, 10, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'visible', 'off');
    hold on;
    p1 = plot(0.1 * (1:length(lbtekData)), lbtekData, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 0.5, 'DisplayName', 'LBTEK');
    p2 = plot(0.1 * (1:length(exceData)), exceData, 'LineStyle', '-', 'Color', 'b', 'LineWidth', 0.5, 'DisplayName', 'Excelitas');
    hold off;

    xlabel('时间 (微秒)');
    ylabel('信号 (光子数)');
    title('1064nm累积58400脉冲后暗计数信号对比');

    xlim([0, 200]);
    ylim([0, 10]);

    text(0.1, 0.8, sprintf('LBTEK暗计数共%d(%dHz);\nExcelitas暗计数共%d(%dHz);\nLBTEK=%4.2fExcelitas', sumDarkCountLBTEK, round(sumDarkCountLBTEK / (100e-9 * length(lbtekData) * 58400)), sumDarkCountExcelitas, round(sumDarkCountExcelitas / (100e-9 * length(exceData) * 58400)), sumDarkCountLBTEK / sumDarkCountExcelitas), 'Units', 'normalized', 'FontSize', 12);

    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'layer', 'top', 'TickLen', [0.02, 0.01], 'FontSize', 12);

    legend([p1, p2], 'Location', 'NorthEast');

    export_fig(gcf, fullfile(outputFolder, sprintf('dark_count_comparison_%02d.png', iFile)), '-r300');
    
    close;
end

%% Efficiency Comparison
subDataPath = fullfile(dataFolder, '\2_连续测试2h_爱尔达');
files = listfile(subDataPath, '.*dat', 1);
allLbtekData = [];
allExceData = [];
mTime = [];

for iFile = 1:length(files)
    thisLidarData = readALADat(files{iFile}, 'nMaxBin', 1900);
    thisLbtekData = thisLidarData.rawSignal(:, 1);
    thisExceData = thisLidarData.rawSignal(:, 2);

    allLbtekData = [allLbtekData, thisLbtekData];
    allExceData = [allExceData, thisExceData];
    mTime = [mTime, thisLidarData.mTime];

    sumDarkCountLBTEK = sum(thisLbtekData);
    sumDarkCountExcelitas = sum(thisExceData);

    figure('Position', [0, 10, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'Visible', 'off');
    hold on;
    p1 = plot(0.1 * (1:length(thisLbtekData)), thisLbtekData / (100e-3*58400), 'LineStyle', '-', 'Color', 'r', 'LineWidth', 0.5, 'DisplayName', 'LBTEK');
    p2 = plot(0.1 * (1:length(thisExceData)), thisExceData / (100e-3*58400), 'LineStyle', '-', 'Color', 'b', 'LineWidth', 0.5, 'DisplayName', 'Excelitas');
    hold off;

    xlabel('时间 (微秒)');
    ylabel('信号 (MHz)');
    title('1064nm累积58400脉冲后大气探测信号对比');

    xlim([0, 200]);
    ylim([1e-4, 10]);

    text(0.1, 0.8, sprintf('LBTEK=%4.2fExcelitas', sum(thisLbtekData(100:500)) / sum(thisExceData(100:500))) , 'Units', 'normalized', 'FontSize', 12);

    set(gca, 'XMinorTick', 'on', 'YTick', 10.^(-4:1), 'YMinorTick', 'on', 'Box', 'on', 'layer', 'top', 'TickLen', [0.02, 0.01], 'FontSize', 12, 'YScale', 'log');

    legend([p1, p2], 'Location', 'NorthEast');

    export_fig(gcf, fullfile(outputFolder, sprintf('signal_comparison_%02d.png', iFile)), '-r300');

    close;
end

%% 展示累计后的信号
meanAllLbtekData = mean(allLbtekData, 2);
meanAllExceData = mean(allExceData, 2);
meanAllLbtekData(meanAllLbtekData <= 0) = NaN;
meanAllExceData(meanAllExceData <= 0) = NaN;
bgAllLbtekData = mean(meanAllLbtekData(1800:1900));
bgAllExceData = mean(meanAllExceData(1800:1900));
snrAllLbtekData = (meanAllLbtekData - bgAllLbtekData) ./ sqrt(meanAllLbtekData);
snrAllExceData = (meanAllExceData - bgAllExceData) ./ sqrt(meanAllExceData);

figure('Position', [0, 10, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'Visible', 'off');
hold on;
p1 = plot(0.1 * (1:length(meanAllLbtekData)), meanAllLbtekData / (100e-3*58400), 'LineStyle', '-', 'Color', 'r', 'LineWidth', 0.5, 'DisplayName', 'LBTEK');
p2 = plot(0.1 * (1:length(meanAllExceData)), meanAllExceData / (100e-3*58400), 'LineStyle', '-', 'Color', 'b', 'LineWidth', 0.5, 'DisplayName', 'Excelitas');
hold off;

xlabel('时间 (微秒)');
ylabel('信号 (MHz)');
title('1064nm累积后大气探测信号对比');

xlim([0, 200]);
ylim([1e-5, 10]);

text(0.1, 0.8, sprintf('探测效率：LBTEK=%4.2fExcelitas;\n背景：LBTEK=%4.2fExcelitas', sum(meanAllLbtekData(100:500)) / sum(meanAllExceData(100:500)), bgAllLbtekData / bgAllExceData), 'Units', 'normalized', 'FontSize', 12);

set(gca, 'XMinorTick', 'on', 'YTick', 10.^(-5:1), 'YMinorTick', 'on', 'Box', 'on', 'layer', 'top', 'TickLen', [0.02, 0.01], 'FontSize', 12, 'YScale', 'log');

legend([p1, p2], 'Location', 'NorthEast');

export_fig(gcf, fullfile(outputFolder, 'all_signal_comparison.png'), '-r300');

close;

figure('Position', [0, 10, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'Visible', 'off');
hold on;
p1 = plot(0.1 * (1:length(meanAllLbtekData)), (meanAllLbtekData - bgAllLbtekData) / (100e-3*58400), 'LineStyle', '-', 'Color', 'r', 'LineWidth', 0.5, 'DisplayName', 'LBTEK');
p2 = plot(0.1 * (1:length(meanAllExceData)), (meanAllExceData - bgAllExceData) / (100e-3*58400), 'LineStyle', '-', 'Color', 'b', 'LineWidth', 0.5, 'DisplayName', 'Excelitas');
hold off;

xlabel('时间 (微秒)');
ylabel('信号 (MHz)');
title('1064nm累积后大气探测信号（扣除背景）对比');

xlim([0, 200]);
ylim([1e-6, 10]);

text(0.1, 0.8, sprintf('LBTEK=%4.2fExcelitas;\n', sum(meanAllLbtekData(100:500)) / sum(meanAllExceData(100:500))), 'Units', 'normalized', 'FontSize', 12);

set(gca, 'XMinorTick', 'on', 'YTick', 10.^(-6:1), 'YMinorTick', 'on', 'Box', 'on', 'layer', 'top', 'TickLen', [0.02, 0.01], 'FontSize', 12, 'YScale', 'log');

legend([p1, p2], 'Location', 'NorthEast');

export_fig(gcf, fullfile(outputFolder, 'all_signal_no_bg_comparison.png'), '-r300');

close;

figure('Position', [0, 10, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'Visible', 'off');
hold on;
p1 = plot(0.1 * (1:length(meanAllLbtekData)), 10 * log10(snrAllLbtekData), 'LineStyle', '-', 'Color', 'r', 'LineWidth', 0.5, 'DisplayName', 'LBTEK');
p2 = plot(0.1 * (1:length(meanAllExceData)), 10 * log10(snrAllExceData), 'LineStyle', '-', 'Color', 'b', 'LineWidth', 0.5, 'DisplayName', 'Excelitas');
hold off;

xlabel('时间 (微秒)');
ylabel('信噪比 (dB)');
title('1064nm累积后大气探测信号信噪比对比');

xlim([0, 200]);
ylim([-20, 20]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'layer', 'top', 'TickLen', [0.02, 0.01], 'FontSize', 12);

legend([p1, p2], 'Location', 'NorthEast');

export_fig(gcf, fullfile(outputFolder, 'all_snr_comparison.png'), '-r300');

close;

%% 展示时空高度图
figure('Position', [0, 10, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'Visible', 'on');

height = (1:size(allLbtekData, 1)) * 15 - 500;

hold on;
p1 = pcolor(mTime, height * 1e-3, allLbtekData .* repmat((height' * 1e-3).^2, 1, size(allLbtekData, 2)) / (100e-3*58400));
p1.EdgeColor = 'none';
hold off;

colormap('jet');
caxis([0, 0.1]);

xlabel('时间（北京时）');
ylabel('高度（千米）');
title('1064nm大气探测时空图 (LBTEK)');

xlim([mTime(1), mTime(end)]);
ylim([0, 10]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickDir', 'out', 'layer', 'top', 'TickLen', [0.02, 0.01], 'FontSize', 12);
datetick('x', 'HH:MM', 'keepticks');

cb = colorbar;
cb.Label.String = '归一化信号 (MHz·km^2)';

export_fig(gcf, fullfile(outputFolder, 'colorplot_lbtek.png'), '-r300');

figure('Position', [0, 10, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'Visible', 'on');

hold on;
p1 = pcolor(mTime, height * 1e-3, allExceData .* repmat((height' * 1e-3).^2, 1, size(allExceData, 2)) / (100e-3*58400));
p1.EdgeColor = 'none';
hold off;

colormap('jet');
caxis([0, 0.1]);

xlabel('时间（北京时）');
ylabel('高度（千米）');
title('1064nm大气探测时空图 (Excelitas)');

xlim([mTime(1), mTime(end)]);
ylim([0, 10]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickDir', 'out', 'layer', 'top', 'TickLen', [0.02, 0.01], 'FontSize', 12);
datetick('x', 'HH:MM', 'keepticks');

cb = colorbar;
cb.Label.String = '归一化信号 (MHz·km^2)';

export_fig(gcf, fullfile(outputFolder, 'colorplot_exec.png'), '-r300');