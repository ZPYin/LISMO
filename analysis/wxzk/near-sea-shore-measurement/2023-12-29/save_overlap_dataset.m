clc;

%% Parameter Definition
overlapFile = 'overlap_20231229.mat';
dataPath = 'E:\wxzk_fog_measurements\RawData\QingDao1\2023\12\29';
savePath = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\wxzk\near-sea-shore-measurement\2023-12-29';

%% Load Data

% overlap
olRes = load(overlapFile);

% lidar
files = listfile(dataPath, '\w*', 2);
liData = readVIS(files{1});
range = ((1:liData.nBins(1)) - 0.5) * liData.hRes(1) - 48.75;

%% Signal Preprocess
bg = mean(squeeze(liData.rawSignal(:, 1, (end - 20):end)), 1);
liSig = squeeze(liData.rawSignal(:, 1, :)) - repmat(bg, liData.nBins(1), 1);

%% Save Data
saveFile = fullfile(savePath, 'overlap_cor_demo_file.txt');
fid = fopen(saveFile, 'w');

fprintf(fid, 'Dataset for the demonstration of overlap correction process.\n');
fprintf(fid, 'Height (m); Overlap factor; Signal\n');

for iH = 1:size(liData.rawSignal, 3)
    fprintf(fid, '%f %f %f\n', range(iH), olRes.ov(iH), liSig(iH));
end

fclose(fid);