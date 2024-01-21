% fog lidar with polarization

clc;
close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\海试数据、效果图\10月21长航实验原始数据';
tRange = [datenum(2023, 10, 21, 17, 0, 0), datenum(2023, 10, 21, 18, 14, 0)];
distOffset = 55.75;
olHeight = 300;
refDist = [1500, 3000];
lr = 30;
maxRange = 10;   % km
fitRange = [0, 4];   % signal glue range. (MHz)
flagReadData = false;
overlapCor = false;
overlapFile = 'overlap_20231021.mat';
visFile = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\wxzk\shipborne-measurements\visibility.mat';

%% Read Data
if flagReadData
    lData = readVIS(dataFolder, 'tRange', tRange, 'isDir', true, 'debug', true);
end

% read visibility
visData = load(visFile);
load('vis_colormap.mat');

%% Preprocessing
range = ((1:lData.nBins(1)) - 0.5) * lData.hRes(1) - distOffset;
sigFR = squeeze(lData.rawSignal(:, 3, :));
sigNR = squeeze(lData.rawSignal(:, 4, :));
sigFRNoBG = sigFR - nanmean(sigFR(:, (end - 15):end), 2);
sigNRNoBG = sigNR - nanmean(sigNR(:, (end - 15):end), 2);
% isSaturated = (sigFR / (lData.hRes(1) / 15 * lData.nShots(1) * 100 * 1e-3) > fitRange(2));
isSaturated = false(size(sigFR));
signal = sigFR;
signal(isSaturated) = ((sigNR(isSaturated) / (lData.hRes(1) / 15 * lData.nShots(1) * 100 * 1e-3)) * 52.9573 - 0.0060) * (lData.hRes(1) / 15 *lData.nShots(1) * 100 * 1e-3);
bg = nanmean(signal(:, (end - 10):end), 2);
signal = signal - repmat(bg, 1, lData.nBins(1));
rcs = signal .* repmat(range, length(lData.hRes), 1).^2;
snr = signal ./ sqrt(signal + repmat(bg, 1, lData.nBins(1)));
if overlapCor
    ol = load(overlapFile);
    rcs = rcs ./ repmat(transpose(ol.ov), size(rcs, 1), 1);
end

%% Xian method
iPrf = 70;
ext = extRet_Xian(range, rcs(iPrf, :) ./ range.^2, bg(iPrf), 'minSNR', 1, 'rangeFullOverlap', 700);

figure;

subplot(311);
plot(range, rcs(iPrf, :));
xlim([0, 5000]);

ylabel('rcs');

subplot(312);
plot(range, ext * 1e3);
xlim([0, 5000]);
ylabel('extinction (km-1)');

subplot(313);
plot(range, ext2vis(ext) * 1e3);
xlim([0, 5000]);
ylabel('visibility (km-1)');
