% This script converts the lidar data from the *.txt/*.dat files to .mat files
%
% Author: Zhenping Yin
% Email: zp.yin@whu.edu.cn
% Date: 2025-01-15

clc;
close all;

%% Parameter Definition
l0Folder = 'G:\backup\vis-lidar\20260428-tianjin-evaluation\前散仪对比数据\VIS1\L0';   % 雷达原始数据文件目录
distOffset = -21;   % 预触发点个数（可通过信号第一个峰值进行判断）
                    % 一般可以先运行代码，读取完数据后，显示其中一根廓线，然后看前面
                    % 的独立峰值信号位置
                    % figure; plot(lidarSig(:, 1));

%% Find Subfolders
L0SubFolders = listdir(l0Folder, '\w*', 1);

% Level 0
for iFolder = 1:length(L0SubFolders)

    fprintf('Finished %6.2f%%: processing %s\n', ...
        (iFolder - 1) / length(L0SubFolders) * 100, L0SubFolders{iFolder});

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

fprintf('Finished 100%%\n');