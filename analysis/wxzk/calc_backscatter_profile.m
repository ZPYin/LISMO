%------------------------------------------------------------------------------%
% Test Jinghong Xian's method for visibility calculation
% 2023-05-05
% Zhenping Yin
% zp.yin@whu.edu.cn
%------------------------------------------------------------------------------%

%% Parameter Definition
dataFolder = 'D:\Data\wxzk_fog_measurements\H_SCAN\2022\8\18\H_SCAN60_270_3_20220818134844';
debug = true;
bgIdx = [2800, 2850];
minRangeOL = 350;

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
data.rcs = transpose(squeeze(data.signal(:, 1, :))) .* repmat(data.distance, size(data.signal, 1), 1).^2;

%% Rayleigh scattering
[temperature, pressure, ~, ~] = read_meteordata(data.startTime, data.height + 0, ...
    'meteor_data', 'standard_atmosphere', ...
    'station', 'xxx');
[mBsc, mExt] = rayleigh_scattering(1064, pressure, temperature, 380, 60);

%% Fernald method

%% Xian
iPrf = 1;

[~, idxMinOL] = min(abs(data.distance - minRangeOL));
SNR = lidarSNR(squeeze(data.signal(iPrf, 1, :)), bg(iPrf, 1), 'bgBins', bgIdx);
RCS = data.rcs(iPrf, :);
RCS_A = RCS(idxMinOL);
ratioAB = RCS / RCS_A;
isInvalid = (data.distance < minRangeOL) | (SNR' <= 0.5);
ratioAB(isInvalid | (ratioAB <= 0)) = NaN;
[~, idxB] = min(ratioAB);
%idxB = 2000;
RCS_B = RCS(idxB);

aExt = NaN(size(RCS));
if (idxB ~= idxMinOL)
    I_A_B = trapz(data.distance(idxMinOL:idxB), RCS(idxMinOL:idxB));
else
    I_A_B = 0;
end
for iBin = (idxMinOL + 1):idxB
    aExt(iBin) = RCS(iBin) / (I_A_B / (1 - RCS_B / RCS_A) - trapz(data.distance(idxMinOL:iBin), RCS(idxMinOL:iBin)));
end

vis = ext2vis(aExt);

figure('Position', [0, 30, 700, 400], 'Units', 'Pixels', 'Color', 'w');
subpos = subfigPos([0.1, 0.12, 0.87, 0.83], 1, 4, 0.03, 0);

subplot('Position', subpos(1, :), 'Units', 'Normalized');
semilogx(RCS, data.distance);
ylim([0, 20000]);
xlim([1e7, 1e11]);

xlabel('R.C.S');
ylabel('Range (m)');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

subplot('Position', subpos(2, :), 'Units', 'Normalized');
plot(SNR, data.distance); hold on;
plot([0.5, 0.5], [0, 20000], '--k');
scatter(SNR(idxB), data.distance(idxB), 'markeredgecolor', 'r');
ylim([0, 20000]);
xlim([1e-1, 1e3]);

xlabel('SNR');
ylabel('');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'XScale', 'log');

subplot('Position', subpos(3, :), 'Units', 'Normalized');
semilogx(aExt, data.distance);
ylim([0, 20000]);
xlim([1e-6, 1e-2]);

xlabel('ext. (m^-1)');
ylabel('');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '');

subplot('Position', subpos(4, :), 'Units', 'Normalized');
plot(vis, data.distance); hold on;
plot([20000, 20000], [0, 20000], '--k');
plot([20000, 20000], [0, 20000], '--k');
ylim([0, 20000]);
xlim([1e2, 1e5]);

xlabel('vis. (m)');
ylabel('');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'xScale', 'log');