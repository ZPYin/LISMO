clc;
%close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\12��27����ͧ��������';
distOffset = -52.5;
olFile = 'overlap_1226.mat';
linFitRange = [1000, 1700];

%% List Data Files
dataFiles = listfile(dataFolder, '\w*.VIS', 1);
dataFiles = dataFiles(177:198);
data = readVIS(dataFiles);

%% Overlap Correction
range = ((1:data.nBins(1)) - 0.5) * data.hRes(1) + distOffset;
bg = nanmean(squeeze(data.rawSignal(:, 3, 2900:2980)), 2);
signal = squeeze(data.rawSignal(:, 3, :)) - repmat(bg, 1, data.nBins(1));
rcs = nansum(signal .* repmat(range, length(data.hRes), 1).^2, 1);
snr = (signal) ./ sqrt(squeeze(data.rawSignal(:, 3, :)));

%% Calculate RCS Slope
isInFit = (range >= linFitRange(1)) & (range <= linFitRange(2));
fitRange = range(isInFit);
fitLogRCS = log(rcs(isInFit));
[offset, slope] = chi2fit(fitRange, fitLogRCS, log(sqrt(rcs(isInFit))));
positiveRCS = rcs;
positiveRCS(positiveRCS <= 0) = NaN;
fullRCS = exp(offset + slope * range);

%% Display
figure('Position', [0, 0, 700, 600], 'Units', 'Pixels', 'Color', 'w');

subplot(211);
p1 = semilogy(range, rcs, 'LineWidth', 2, 'DisplayName', 'ԭʼ����');
hold on;
p2 = semilogy(range, fullRCS, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', '�ؽ���');

xlabel('���� (��)');
ylabel('���������ź�');
title('2023��12��26�� 12:00');

xlim([0, 2000]);
ylim([5e8, 1e12]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
legend([p1, p2], 'Location', 'NorthEast');

subplot(212);
p1 = plot(range, rcs ./ fullRCS, '-k', 'LineWidth', 2);
hold on;
p2 = plot([0, 1e5], [1, 1], '--r');

xlabel('���� (��)');
ylabel('�ص�����');

xlim([0, 2000]);
ylim([-0.1, 1.2]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

%% Save Overlap
ov = rcs ./ fullRCS;
ov(range >= 600) = 1;
save(olFile, 'ov', 'range');