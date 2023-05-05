%------------------------------------------------------------------------------%
% Test cloud detection algorithm
% 2023-05-05
% Zhenping Yin
% zp.yin@whu.edu.cn
%------------------------------------------------------------------------------%

global LISMO_VARS;

%% Parameter Definition
dataFolder = 'D:\Data\wxzk_fog_measurements\H_SCAN\2022\8\18\H_SCAN60_270_3_20220818225724';
debug = true;
bgIdx = [2800, 2900];
tRes = datenum(0, 1, 0, 0, 1, 0);

%% Read Data
dataFiles = listfile(dataFolder, '.*.VIS', 2);
data = readVIS(dataFiles, 'debug', debug);

%% Signal Preprocess
data.height = (1:data.nBins) * data.hRes(1) * cos(data.zenithAng(1) / 180 * pi);
data.distance = (1:data.nBins) * data.hRes(1);

% remove background
isBgBins = false(1, length(data.height));
isBgBins(bgIdx(1):bgIdx(2)) = true;
bg = nanmean(data.rawSignal(:, :, isBgBins), 3);
data.signal = data.rawSignal - repmat(bg, 1, 1, size(data.rawSignal, 3));

% range corrected signal
data.rcs = squeeze(data.signal(:, 1, :)) .* repmat(data.distance, size(data.signal, 1), 1).^2;

%% cloud detection with aerosol
mTime = data.startTime;
iPrf = 3;
minDepth = data.hRes(iPrf) * cos(data.zenithAng(iPrf) / 180 * pi) * 3;
[~, layerStatus] = cloudDetect_Zhao(mTime(iPrf), data.height, squeeze(data.signal(iPrf, 1, :)), bg(iPrf, 1), 'minDepth', minDepth, 'detectRange', [0, 1500], 'smoothWin', 4, 'minSNR', 3, 'heightFullOverlap', 300 * cos(data.zenithAng(iPrf) / 180 * pi));
rcs = transpose(squeeze(data.signal(iPrf, 1, :))) .* data.distance.^2;

% show feature detection results of single profile
figure('Position', [0, 30, 500, 500], 'Units', 'Pixels', 'Color', 'w');

subPos = subfigPos([0.1, 0.1, 0.87, 0.87], 1, 2, 0.05, 0);

subplot('Position', subPos(1, :), 'Units', 'Normalized', 'Color', 'w');
semilogx(rcs, data.distance / 1e3); hold on;
plot(rcs(layerStatus == 1), data.distance(layerStatus == 1) / 1e3, 'Color', 'r', 'LineStyle', 'none', 'Marker', '.');
plot(rcs(layerStatus == 2), data.distance(layerStatus == 2) / 1e3, 'Color', 'g', 'LineStyle', 'none', 'Marker', '.');

xlim([1e5, 1e12]);
ylim([0, 22]);

xlabel('range-cor. sig. [a.u]');
ylabel('distance (km)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

subplot('Position', subPos(2, :), 'Units', 'Normalized', 'Color', 'w');

plot(layerStatus, data.distance / 1e3); hold on;

xlim([-0.5, 2.5]);
ylim([0, 22]);

xlabel('feature');
ylabel('');

set(gca, 'XTick', 0:1:2, 'XTickLabel', {'none', 'cloud', 'aerosol'}, 'YTicklabel', '', 'Box', 'on', 'YMinorTick', 'on');

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'results', 'wxzk', 'cloud_detection_with_aerosol.png'), '-r300');

%% cloud detection without aerosol
mTime = data.startTime;
iPrf = 1;
minDepth = data.hRes(iPrf) * cos(data.zenithAng(iPrf) / 180 * pi) * 3;
[~, layerStatus] = cloudDetect_Zhao(mTime(iPrf), data.height, squeeze(data.signal(iPrf, 1, :)), bg(iPrf, 1), 'minDepth', minDepth, 'detectRange', [0, 1500], 'smoothWin', 4, 'minSNR', 3, 'heightFullOverlap', 300 * cos(data.zenithAng(iPrf) / 180 * pi));
rcs = transpose(squeeze(data.signal(iPrf, 1, :))) .* data.distance.^2;

% show feature detection results of single profile
figure('Position', [0, 30, 500, 500], 'Units', 'Pixels', 'Color', 'w');

subPos = subfigPos([0.1, 0.1, 0.87, 0.87], 1, 2, 0.05, 0);

subplot('Position', subPos(1, :), 'Units', 'Normalized', 'Color', 'w');
semilogx(rcs, data.distance / 1e3); hold on;
plot(rcs(layerStatus == 1), data.distance(layerStatus == 1) / 1e3, 'Color', 'r', 'LineStyle', 'none', 'Marker', '.');
plot(rcs(layerStatus == 2), data.distance(layerStatus == 2) / 1e3, 'Color', 'g', 'LineStyle', 'none', 'Marker', '.');

xlim([1e5, 1e12]);
ylim([0, 22]);

xlabel('range-cor. sig. [a.u]');
ylabel('distance (km)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

subplot('Position', subPos(2, :), 'Units', 'Normalized', 'Color', 'w');

plot(layerStatus, data.distance / 1e3); hold on;

xlim([-0.5, 2.5]);
ylim([0, 22]);

xlabel('feature');
ylabel('');

set(gca, 'XTick', 0:1:2, 'XTickLabel', {'none', 'cloud', 'aerosol'}, 'YTicklabel', '', 'Box', 'on', 'YMinorTick', 'on');

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'results', 'wxzk', 'cloud_detection.png'), '-r300');

%% cloud detection without cloud
mTime = data.startTime;
iPrf = 4;
minDepth = data.hRes(iPrf) * cos(data.zenithAng(iPrf) / 180 * pi) * 3;
[~, layerStatus] = cloudDetect_Zhao(mTime(iPrf), data.height, squeeze(data.signal(iPrf, 1, :)), bg(iPrf, 1), 'minDepth', minDepth, 'detectRange', [0, 1500], 'smoothWin', 4, 'minSNR', 3, 'heightFullOverlap', 300 * cos(data.zenithAng(iPrf) / 180 * pi));
rcs = transpose(squeeze(data.signal(iPrf, 1, :))) .* data.distance.^2;

% show feature detection results of single profile
figure('Position', [0, 30, 500, 500], 'Units', 'Pixels', 'Color', 'w');

subPos = subfigPos([0.1, 0.1, 0.87, 0.87], 1, 2, 0.05, 0);

subplot('Position', subPos(1, :), 'Units', 'Normalized', 'Color', 'w');
semilogx(rcs, data.distance / 1e3); hold on;
plot(rcs(layerStatus == 1), data.distance(layerStatus == 1) / 1e3, 'Color', 'r', 'LineStyle', 'none', 'Marker', '.');
plot(rcs(layerStatus == 2), data.distance(layerStatus == 2) / 1e3, 'Color', 'g', 'LineStyle', 'none', 'Marker', '.');

xlim([1e5, 1e12]);
ylim([0, 22]);

xlabel('range-cor. sig. [a.u]');
ylabel('distance (km)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

subplot('Position', subPos(2, :), 'Units', 'Normalized', 'Color', 'w');

plot(layerStatus, data.distance / 1e3); hold on;

xlim([-0.5, 2.5]);
ylim([0, 22]);

xlabel('feature');
ylabel('');

set(gca, 'XTick', 0:1:2, 'XTickLabel', {'none', 'cloud', 'aerosol'}, 'YTicklabel', '', 'Box', 'on', 'YMinorTick', 'on');

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'results', 'wxzk', 'cloud_detection_without_cloud.png'), '-r300');