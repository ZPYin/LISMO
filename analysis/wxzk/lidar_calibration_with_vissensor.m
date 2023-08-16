clc;

%% Parameter Definition
lidarDataPath = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\RawData\QingDao\2023\8\10\H_SCAN30_150_2_20230810204742';
visFilePath = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\visibility.txt';
tRange = [datenum(2023, 8, 10, 20, 47, 0), datenum(2023, 8, 10, 20, 57, 0)];
overlapCor = false;

%% Read Data

% read lidar data
dataFiles = listfile(lidarDataPath, '\w*', 1);
dataFilesInTimeRange = cell(0);
for iFile = 1:length(dataFiles)
    [~, tmp, ext] = fileparts(dataFiles{iFile});
    fileTime = datenum(tmp(3:end), 'yyyymmddHHMMSS');

    if (fileTime >= tRange(1)) && (fileTime <= tRange(2))
        dataFilesInTimeRange = cat(2, dataFilesInTimeRange, dataFiles{iFile});
    end
end

lidarData = readVIS(dataFilesInTimeRange);
range1 = ((1:lidarData.nBins(1)) - 0.5) * lidarData.hRes(1) - 48.75;
bg = nanmean(squeeze(lidarData.rawSignal(:, 1, 2900:2980)), 2);
signal = squeeze(lidarData.rawSignal(:, 1, :)) - repmat(bg, 1, lidarData.nBins(1));
if overlapCor
    load('overlap_0527.mat');
    signal = signal ./ repmat(ov, size(signal, 1), 1);
end
rcs = signal .* repmat(range1, length(lidarData.hRes), 1).^2;

% read visibility snesor
[~, ~, ext] = fileparts(visFilePath);

switch ext
case '.txt'
    fid = fopen(visFilePath, 'r');

    data = textscan(fid, '%s%f', 'delimiter', ',');

    fclose(fid);
    visTime = [];
    vis = [];
    for iRow = 1:length(data{1})
        visTime = cat(2, visTime, datenum(data{1}{iRow}, 'yyyy-mm-dd HH:MM:SS'));
        vis = cat(2, vis, data{2}(iRow));
    end
case '.xls'
    [~, ~, visData] = xlsread(visFilePath, 'Êý¾Ý');
    visTime = [];
    vis = [];
    for iRow = 2:size(visData, 1)
        tmp = datenum(['2023', visData{iRow, 2}(1:2), visData{iRow, 2}(4:5), visData{iRow, 2}(7:8), visData{iRow, 2}(10:11)], 'yyyymmddHHMM');
    
        if (tmp >= tRange(1)) && (tmp <= tRange(2))
            visTime = cat(2, visTime, tmp);
            vis = cat(2, vis, str2num(visData{iRow, 4}));
        end
    end
otherwise
    error('Unknown file type.');
end

%% Display
figure('Position', [0, 0, 500, 280], 'color', 'w', 'visible', 'on');
p1 = pcolor(lidarData.startTime, range / 1e3, transpose(rcs));
p1.EdgeColor = 'none';

colormap('jet');

xlabel('Local Time');
ylabel('Distance (km)');

xlim([lidarData.startTime(1), lidarData.startTime(end)]);
ylim([0, 10]);
caxis([0, 6e8]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickDir', 'out', 'Box', 'on');
datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

colorbar();

figure('color', 'w', 'visible', 'on');
plot(visTime, vis / 1e3, '-k');

xlabel('Local Time');
ylabel('Visibility (km)');

xlim([lidarData.startTime(1), lidarData.startTime(end)]);
ylim([0, 40]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');
datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

%% Calibration
height = zeros(size(range));
[mBsc, mExt] = MolModel(height, 1064, 'meteor', 'standard_atmosphere');

[~, visIdx] = min(abs(visTime - mean(lidarData.startTime)));
refVis = vis(visIdx);
refExt = vis2ext(refVis);

% [aBsc1, ~] = fernald(range, mean(signal, 1), mean(bg), 45, [970, 1030], refExt / 45, mBsc);
[aBsc1, ~] = fernald(range, mean(signal, 1), mean(bg), 45, [9700, 10000], 1.3e-4 / 45, mBsc);


%% Display
figure;

subplot(211);
hold on;
semilogy(range, mean(rcs, 1));
semilogy(range, 2.4e15*mBsc .* exp(-2 * mExt .* [range(1), diff(range)]));
hold off;

xlim([0, 10000]);
set(gca, 'YScale', 'log');

subplot(212);
hold on;
semilogy(range, aBsc1.*45);

%ylim([]);
xlim([0, 10000]);

set(gca, 'YScale', 'linear');

%% Output
lc = mean(rcs) / (refExt / 45);
figure;
hold on;
plot(range, lc);

title(datestr(mean(lidarData.startTime), 'yyyy-mm-dd HH:MM'));