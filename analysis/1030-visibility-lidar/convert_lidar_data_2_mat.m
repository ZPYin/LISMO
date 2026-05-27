% This script converts the lidar data from the .txt files to .mat files
% Author: Zhenping Yin
% Date: 2025-01-15

clc;
close all;

%% Parameter Definition
l0Folder = 'G:\backup\vis-lidar\20251221-wuhan\L0';   % 雷达原始数据文件目录
l1Folder = 'G:\backup\vis-lidar\20251221-wuhan\L1';   % 雷达产品文件目录
distOffset = -21;   % 预触发点个数（可通过信号第一个峰值进行判断）

%% Find Subfolders
L0SubFolders = listdir(l0Folder, '\w*', 1);
L1SubFolders = listdir(l1Folder, '\w*', 1);

% Level 0
for iFolder = 1:length(L0SubFolders)

    fprintf('Finished %6.2f%%: processing %s\n', (iFolder - 1) / length(L0SubFolders) * 100, L0SubFolders{iFolder});

    % Find Files
    thisDate = datenum(basename(L0SubFolders{iFolder}), 'yyyymmdd');
    lData1 = readALADats(L0SubFolders{iFolder}, 'datatype', 2);

    lidarSig = squeeze(lData1.rawSignal);
    range = transpose(((1:size(lData1.rawSignal, 2)) - 0.5 + distOffset) * lData1.hRes(1));
    mTime = lData1.mTime;

    % Save to .mat file
    save(fullfile(l0Folder, sprintf('%s_vis_lidar_l0.mat', datestr(thisDate, 'yyyy-mm-dd'))), ...
        'lidarSig', 'range', 'mTime');

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