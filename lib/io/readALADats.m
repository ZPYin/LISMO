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
%    datatype: numeric
%        1: origianl output from ALA data recorder
%        2: renaming to fulfill CMA standard
%        (default: 1)
%    flagSort: logical
%        whether to sort data in time.
%
% OUTPUTS:
%    data: struct
%        rawSignal: matrix (channel x height x time)
%        mTime: numeric
%        nPretrigger: numeric
%        hRes: numeric
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
addParameter(p, 'flagSort', false, @islogical);
addParameter(p, 'nMaxBin', [], @isnumeric);
addParameter(p, 'debug', false, @islogical);
addParameter(p, 'datatype', 1, @isnumeric);

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
        fprintf('Finish %6.2f%%\n', (iFile - 1) / length(dataFiles) * 100);
    end

    thisFilename = basename(dataFiles{iFile});
    switch p.Results.datatype
        case 1
            mTime = datenum(thisFilename((end-16):(end-4)), 'yymmdd-HHMMSS');
        case 2
            mTime = datenum(thisFilename(16:29), 'yyyymmddHHMMSS');
        otherwise
            error('Unsupported timestamp type: %d', p.Results.timestamp_type);
    end
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

if p.Results.flagSort
    [~, idxS] = sort(data.mTime);

    data.mTime = data.mTime(idxS);
    data.rawSignal = data.rawSignal(:, :, idxS);
    data.channelLabel = data.channelLabel;
    data.hRes = thisData.hRes;
    data.nShots = data.nShots(:, idxS);
end

end