% Convert GBQ L2 data to mat file
% Author: Zhenping Yin
% Date: 2025-03-13

clc;
close all;

%% Parameter Definition
L2Folder = 'C:\Users\ZPYin\Documents\data\CMA-Data-Backup\I_54511';
saveFolder = 'C:\Users\ZPYin\Documents\data\CMA-Data-Backup\I_54511';
nBin = 4000;

%% Find Subfolders
L2SubFolders = listdir(L2Folder, '\w{8}', 1);
for iFolder = 1:length(L2SubFolders)

    % Find Files
    thisDate = datenum(L2SubFolders{iFolder}, 'yyyymmdd');
    thisFiles = listfile(fullfile(L2SubFolders{iFolder}, 'LEVEL2_PRODUCT'), 'Z_RADA\w*', 1);

    if isempty(thisFiles)
        continue;
    end

    % Read Lidar Data
    range = zeros(1, nBin);
    mTime = zeros(1, length(thisFiles));
    vis_gbq = zeros(nBin, length(thisFiles));
    for iFile = 1:length(thisFiles)
        fprintf('Finished %6.2f%% for %s: reading %s\n', (iFile - 1) / length(thisFiles) * 100, basename(L2SubFolders{iFolder}), basename(thisFiles{iFile}));

        fid = fopen(thisFiles{iFile}, 'r');
        thisStr = fgetl(fid);
        mTime(iFile) = datenum(thisStr((end - 18):end), 'yyyy-mm-dd HH:MM:SS');

        data = textscan(fid, '%f%f', 'Delimiter', ',', 'HeaderLines', 5);
        fclose(fid);

        range = data{1};
        vis_gbq(:, iFile) = data{2};
    end

    % Save to .mat file
    save(fullfile(saveFolder, sprintf('%s_vis_lidar_l2.mat', datestr(mTime(1), 'yyyy-mm-dd'))), ...
        'range', 'mTime', 'vis_gbq');
end