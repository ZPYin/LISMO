function [data] = readALADats(folder, varargin)
% READALADATS read ALA dat files from a folder.
%
% USAGE:
%    [data] = readALADats(folder)
%
% INPUTS:
%    folder: char
%        data folder.
%
% KEYWORDS:
%    tRange: numeric
%        temporal range for data files.
%    nMaxBin: numeric
%        maximum data bins to load in.
%
% OUTPUTS:
%    data: struct
%        rawSignal: matrix (channel x height x time)
%        mTime: numeric
%        nPretrigger: numeric
%        nShots: numeric
%        channelLabel: cell
%
% HISTORY:
%    2023-07-19: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'folder', @ischar);
addParameter(p, 'tRange', [], @isnumeric);
addParameter(p, 'nMaxBin', [], @isnumeric);
addParameter(p, 'debug', false, @islogical);

parse(p, folder, varargin{:});

%% search data files
dataFiles = listfile(folder, '.*dat', 1);

%% read data
data = struct();
data.mTime = [];
data.rawSignal = [];
data.channelLabel = cell(0);
data.hRes = [];
data.nPretrigger = [];
data.nShots = [];

for iFile = 1:length(dataFiles)
    if p.Results.debug
        fprintf('Finish %6.2f%% ', (iFile - 1) / length(dataFiles) * 100);
    end

    mTime = datenum(dataFiles{iFile}((end-16):(end-4)), 'yymmdd-HHMMSS');
    if ~ isempty(p.Results.tRange)
        isNotInSearchTRange = (mTime < p.Results.tRange(1)) || (mTime > p.Results.tRange(2));

        if isNotInSearchTRange
            continue;
        end
    end

    thisData = readALADat(dataFiles{iFile}, varargin{:});

    data.mTime = cat(2, data.mTime, mTime);
    data.rawSignal = cat(3, data.rawSignal, reshape(transpose(thisData.rawSignal), size(thisData.rawSignal, 2), size(thisData.rawSignal, 1), 1));
    data.channelLabel = thisData.channelLabel;
    data.hRes = thisData.hRes;
    data.nShots = cat(2, data.nShots, thisData.nShots);
end

end