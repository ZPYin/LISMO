% 对比评估麓邦光电和Excelitas公司单光子探测器的性能
% 2025-09-16

%% Parameter Definition
dataFolder = 'C:\Users\zhenp\Documents\Industry\麓邦光电\LBTEK_APD_new_2026\5min-62um-fiber';
outputFolder = 'C:\Users\zhenp\Documents\Industry\麓邦光电\LBTEK_APD_new_2026';

%% Efficiency Comparison
files = listfile(dataFolder, '.*dat', 1);
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
