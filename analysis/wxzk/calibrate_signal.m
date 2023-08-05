clc;
close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\RawData\QingDao\2023\7\1\H_SCAN30_150_2_20230701200244';
distOffset = 48.75;
olHeight = 300;
refH = [10000, 10500];
ratio = 7.5e16;
lr = 50;

%% List Data Files
dataFiles = listfile(dataFolder, '\w*.VIS', 1);
dataFiles = dataFiles(1:35);
data = readVIS(dataFiles);

%% Overlap Correction
range = ((1:data.nBins(1)) - 0.5) * data.hRes(1) - distOffset;
height = zeros(size(range));
bg = nanmean(squeeze(data.rawSignal(:, 1, 2900:2980)), 2);
signal = squeeze(data.rawSignal(:, 1, :)) - repmat(bg, 1, data.nBins(1));
rcs = nanmean(signal .* repmat(range, length(data.hRes), 1).^2, 1);
snr = (signal) ./ sqrt(squeeze(data.rawSignal(:, 1, :)));

%% Rayleigh Scattering
[mBsc, mExt] = MolModel(height, 1064, 'meteor', 'standard_atmosphere');

%% Rayleigh Fit
figure;
semilogy(range, rcs); hold on;
semilogy(range, mBsc .* exp(-2 * cumsum(mExt .* [range(1), diff(range)])) * ratio);

%% Fernald
[aBsc1, ~] = fernald(range, mean(signal, 1), bg, lr, refH, 0, mBsc);
aBsc1(range < olHeight) = aBsc1(find(range > olHeight, 1, 'first'));

figure;
plot(range, aBsc1);

%% Lidar calibration
lc = rcs ./ ((aBsc1 + mBsc) .* exp(-2 * cumsum(-2 * (aBsc1 .* lr + mBsc) .* [range(1), diff(range)])));

figure; 
plot(range, smooth(lc, 15));