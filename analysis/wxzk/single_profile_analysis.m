%------------------------------------------------------------------------------%
% Check signal profile for data quality control.
% 2023-05-05
% Zhenping Yin
% zp.yin@whu.edu.cn
%------------------------------------------------------------------------------%

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\RawData\QingDao\2023\4\19\H_SCAN0_120_2_20230419143839';
debug = true;

%% Read Data
dataFiles = listfile(dataFolder, '.*.VIS', 1);
data = readVIS(dataFiles, 'debug', debug);

%% Data Visualization

%% Single Profile
iPrf = 10;
distArr = (1:data.nBins) * data.hRes(iPrf);
liSig = squeeze(data.rawSignal(iPrf, 1, :));
bgBins = [2800, 2900];
bg = mean(liSig(bgBins(1):bgBins(2)));
liSigNoBg = liSig - bg;
rcs = liSigNoBg .* transpose(distArr).^2;

figure('Name', 'Single Profile', 'Position', [0, 30, 500, 500], 'Units', 'Pixels', 'Color', 'w');

subPos = subfigPos([0.1, 0.1, 0.86, 0.86], 1, 2, 0.04, 0);

subplot('Position', subPos(1, :), 'units', 'Normalized');
semilogx(liSig, distArr, '-k'); hold on;
semilogx(liSigNoBg, distArr, '-r');
xlabel('sig. [a.u.]');
ylabel('disttance (m)');
xlim([1e-1, 1e5]);
title(sprintf('%s', datestr(data.startTime(iPrf), 'yyyy-mm-dd HH:MM:SS')));
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

subplot('Position', subPos(2, :), 'units', 'Normalized');
semilogx(rcs, distArr, '-k');
xlim([5e7, 1e12]);
xlabel('Range Cor. Sig.');
ylabel('');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'Box', 'on');

%% Averaged Profile
liSig = data.rawSignal(iPrf, 1, :);
liSigAvg = mean(data.rawSignal(:, 1, :), 1);
bg = mean(squeeze(data.rawSignal(:, 1, 2800:2900)), 2);
snr = (squeeze(data.rawSignal(:, 1, :)) - repmat(bg, 1, size(data.rawSignal, 3))) ./ sqrt(squeeze(data.rawSignal(:, 1, :)));

figure('Position', [0, 30, 700, 500], 'Units', 'Pixels', 'Color', 'w');
subPos = subfigPos([0.1, 0.1, 0.86, 0.86], 1, 3, 0.04, 0);

subplot('Position', subPos(1, :), 'Units', 'Normalized');
for iPrf = 1:size(data.rawSignal, 1)
    distArr = (1:data.nBins) * data.hRes(iPrf) - 51;
    semilogx(squeeze(data.rawSignal(iPrf, 1, :)), distArr); hold on;
end
xlabel('raw sig. [a.u.]');
ylabel('disttance (m)');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

subplot('Position', subPos(2, :), 'Units', 'Normalized');
semilogx(squeeze(liSigAvg), distArr);
ylabel('');
xlabel('raw sig. [a.u.]');
xlim([1, 1e5]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'Box', 'on');

subplot('Position', subPos(3, :), 'Units', 'Normalized');
semilogx(snr, distArr);
ylabel('');
xlabel('SNR');
xlim([1, 1e5]);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'Box', 'on');