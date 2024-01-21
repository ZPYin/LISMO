clc;
close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\海试数据、效果图\10月21长航实验原始数据';
% tRange = [datenum(2023, 10, 21, 20, 0, 0), datenum(2023, 10, 21, 21, 14, 0)];
tRange = [datenum(2023, 10, 21, 11, 0, 0), datenum(2023, 10, 21, 13, 14, 0)];
fitRange = [0, 4];   % signal linearity range. (MHz)

%% Read Data
data = readVIS(dataFolder, 'tRange', tRange, 'isDir', true, 'debug', true);

%% Data Preprocess
bg = nanmean(data.rawSignal(:, :, (end - 20):(end - 2)), 3);
signal = data.rawSignal - repmat(bg, 1, 1, data.nBins(1));
snr = signal ./ sqrt(data.rawSignal);

frChPC = reshape(squeeze(signal(:, 3, :)), 1, []);
nrChPC = reshape(squeeze(signal(:, 4, :)), 1, []);
frChPoints = frChPC / (data.hRes(1) / 15 * data.nShots(1) * 100 * 1e-3);
nrChPoints = nrChPC / (data.hRes(1) / 15 * data.nShots(1) * 100 * 1e-3);

isInFitRange = (frChPoints >= fitRange(1)) & (frChPoints <= fitRange(2)) & (nrChPoints >= fitRange(1)) & (nrChPoints <= fitRange(2));

[a, b] = chi2fit(nrChPoints(isInFitRange), frChPoints(isInFitRange), sqrt(nrChPC(isInFitRange)) / (data.hRes(1) / 15 * data.nShots(1) * 100 * 1e-3));

fprintf('far-range = near-range * %s + %s\n', b, a);

%% Signal Scatters
figure('Position', [0, 10, 450, 350], 'Units', 'Pixels', 'Color', 'w');

hold on;
s1 = scatter(frChPoints, nrChPoints, 5, 'Marker', '.', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot([0, 10], ([0, 10] - a) / b, '--b', 'LineWidth', 2);
plot([fitRange(2), fitRange(2)], [0, 100], '--m');
plot([0, 100], [fitRange(2), fitRange(2)], '--m');
hold off;

xlabel('Far-Range Signal (MHz)');
ylabel('Near-Range Signal (MHz)');

xlim([0, 6]);
ylim([0, 0.3]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

%% Signal Glue
sigFR = squeeze(signal(1, 3, :)) / (data.hRes(1) / 15 * data.nShots(1) * 100 * 1e-3);
sigNR = squeeze(signal(1, 4, :)) / (data.hRes(1) / 15 * data.nShots(1) * 100 * 1e-3);
isSaturated = (sigFR >= fitRange(2));
sigFRGlue = sigFR;
sigFRGlue(isSaturated) = sigNR(isSaturated) * b + a;

figure('Position', [0, 20, 500, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
p3 = plot(sigFRGlue, '-r', 'Marker', '.', 'MarkerSize', 10, 'DisplayName', 'Glue');
p1 = plot(sigFR, '-k', 'DisplayName', 'Far-Range');
p2 = plot(sigNR, '-b', 'DisplayName', 'Near-Range');

plot([0, 1e5], [fitRange(2), fitRange(2)], '--m');
hold off;

xlabel('Range Bin');
ylabel('Signal (MHz)');

xlim([0, 150]);
ylim([1e-2, 150]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTick', 10.^(-2:1:2), 'Box', 'on', 'YScale', 'log', 'fontsize', 11);

legend([p1, p2, p3], 'Location', 'NorthEast');