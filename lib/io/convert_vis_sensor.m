function [ dataFile ] = convert_vis_sensor(dataOriPath, dataDstPath, varargin)
% CONVERT_VIS_SENSOR convert visibility sensor data to mat.
%
% USAGE:
%    [dataFile] = convert_vis_sensor(dataOriPath, dataDstPath)
%
% INPUTS:
%    dataOriPath: char
%        original data path.
%    dataDstPath: char
%        destination data path.
%
% OUTPUTS:
%    dataFile: char
%
% HISTORY:
%    2023-11-20: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'dataOriPath', @ischar);
addRequired(p, 'dataDstPath', @ischar);
addParameter(p, 'debug', false, @islogical);

parse(p, dataOriPath, dataDstPath, varargin{:});

files = listfile(dataOriPath, '\w*.csv', 1);

vis35 = struct();
vis35.mTime = [];
vis35.vis1min = [];
vis35.vis10min = [];
vis35.temperature = [];
for iFile = 1:length(files)

    if p.Results.debug
        fprintf('Reading %s\n', files{iFile});
    end

    data = readV35(files{iFile});

    vis35.mTime = cat(2, vis35.mTime, data.mTime);
    vis35.vis1min = cat(2, vis35.vis1min, data.vis1min);
    vis35.vis10min = cat(2, vis35.vis10min, data.vis10min);
    vis35.temperature = cat(2, vis35.temperature, data.temp);
end

dataFile = fullfile(dataDstPath, 'vis35.mat');
save(dataFile, 'vis35');

end