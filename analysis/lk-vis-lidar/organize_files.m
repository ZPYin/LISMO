%% Parameter Definition
dataPath = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs';
savePath = 'C:\Users\ZPYin\Documents\Data\CMA-vis-lidar-assessment\highway-obs\results';

%% List Files
dataFiles = listfile(dataPath, '\w*pb', 1);

%% Parse File Info
fileTimeArr = [];
for iFile = 1:length(dataFiles)
    thisFilename = basename(dataFiles{iFile});
    thisFileTime = datenum(thisFilename(1:(end - 3)), 'yyyymmddHHMMSS');

    fileTimeArr = cat(2, fileTimeArr, thisFileTime);
end

%% Create Folders
for iFile = 1:length(dataFiles)
    subFolder = fullfile(savePath, datestr(fileTimeArr(iFile), 'yyyy-mm-dd'));
    if ~ exist(subFolder, 'dir')
        warning('Create folder for saving: %s', subFolder);
        mkdir(subFolder);
    end

    movefile(dataFiles{iFile}, subFolder);
end

fprintf('Finish!\n');