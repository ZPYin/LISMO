% This script converts the lidar data from the .txt files to .mat files
% Author: Zhenping Yin
% Date: 2025-01-15

clc;
close all;

%% Parameter Definition
l0Folder = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\L0';   % 雷达原始数据文件目录
l1Folder = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\L1';   % 雷达产品文件目录
nPretrigger = 12;

%% Find Subfolders
L0SubFolders = listdir(l0Folder, '\w*', 1);
L1SubFolders = listdir(l1Folder, '\w*', 1);

% Level 0
for iFolder = 1:length(L0SubFolders)

    fprintf('Finished %6.2f%%: processing %s\n', (iFolder - 1) / length(L0SubFolders) * 100, L0SubFolders{iFolder});

    % Find Files
    thisDate = datenum(basename(L0SubFolders{iFolder}), 'yyyy-mm-dd');
    thisFiles = listfile(L0SubFolders{iFolder}, '\w*pb', 1);

    if isempty(thisFiles)
        continue;
    end

    % Read Lidar Data
    lidarSig = [];
    laserTemp = [];
    mTime = [];
    for iFile = 1:length(thisFiles)
        % fprintf('Finished %6.2f%%: reading %s\n', (iFile - 1) / length(thisFiles) * 100, thisFiles{iFile});
        thisData = readVisLidarL0(thisFiles{iFile});
        thisData.sig = thisData.sig(nPretrigger:end);

        lidarSig = cat(2, lidarSig, thisData.sig);
        mTime = cat(2, mTime, thisData.startTime);
        laserTemp = cat(2, laserTemp, thisData.laserTemp);
    end

    % Save to .mat file
    save(fullfile(l0Folder, sprintf('%s_vis_lidar_l0.mat', datestr(thisDate, 'yyyy-mm-dd'))), ...
        'lidarSig', 'laserTemp', 'mTime');

end

% Level 1
for iFolder = 1:length(L1SubFolders)

    fprintf('Finished %6.2f%%: processing %s\n', (iFolder - 1) / length(L1SubFolders) * 100, L1SubFolders{iFolder});

    % Find Files
    thisDate = datenum(basename(L1SubFolders{iFolder}), 'yyyy-mm-dd');
    thisFiles = listfile(L1SubFolders{iFolder}, '\w*pb', 1);

    if isempty(thisFiles)
        continue;
    end

    % Read Lidar Data
    ext_LK = [];
    vis_LK = [];
    height_LK = [];
    mTime_LK = [];
    for iFile = 1:length(thisFiles)
        % fprintf('Finished %6.2f%%: reading %s\n', (iFile - 1) / length(thisFiles) * 100, thisFiles{iFile});

        thisL1 = readVisLidarL1(thisFiles{iFile});
        vis_LK = cat(2, vis_LK, thisL1.vis);
        ext_LK = cat(2, ext_LK, thisL1.extinction);
        mTime_LK = cat(2, mTime_LK, thisL1.mTime);
        height_LK = thisL1.height;
    end

    % Save to .mat file
    save(fullfile(l1Folder, sprintf('%s_vis_lidar_l1.mat', datestr(thisDate, 'yyyy-mm-dd'))), ...
        'ext_LK', 'vis_LK', 'height_LK', 'mTime_LK');

end

fprintf('Finished 100%%\n');