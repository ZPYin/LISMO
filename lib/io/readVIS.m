function [data] = readVIS(filename, varargin)
% READVIS Read fog lidar from WXZK.
%
% USAGE:
%    [data] = readVIS(filename)
%
% INPUTS:
%    filename: char | cell
%        absolute path of data file(s).
%
% KEYWORDS:
%    debug: logical
%        whether to start Debug mode. (default: false)
%
% OUTPUTS:
%    data: struct
%        hRes
%        nShots
%        nBins
%        startTime
%        stopTime
%        rawSignal
%        isRain
%        temperature
%        azimuthAng
%        zenithAng
%        latitude
%        longitude
%        asl
%        speed
%
% HISTORY:
%    2023-02-07: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'filename', @(x) ischar(x) || iscell(x));
addParameter(p, 'debug', false, @islogical);

parse(p, filename, varargin{:});

if ischar(filename)
    files = {filename};
else
    files = filename;
end

data = struct();
data.hRes = [];
data.nShots = [];
data.nBins = [];
data.startTime = [];
data.stopTime = [];
data.rawSignal = [];
data.isRain = [];
data.temperature = [];
data.azimuthAng = [];
data.zenithAng = [];
data.latitude = [];
data.longitude = [];
data.asl = [];
data.speed = [];

for iFile = 1:length(files)
    if (p.Results.debug)
        fprintf('Finished %6.5f%%: reading %s\n', (iFile - 1) / length(files) * 100, files{iFile});
    end

    singlePrfData = readSingleVIS(files{iFile}, varargin{:});

    data.hRes = cat(2, data.hRes, singlePrfData.hRes);
    data.nShots = cat(2, data.nShots, singlePrfData.nShots);
    data.nBins = cat(2, data.nBins, singlePrfData.nBins);
    data.startTime = cat(2, data.startTime, singlePrfData.startTime);
    data.stopTime = cat(2, data.stopTime, singlePrfData.stopTime);
    data.rawSignal = cat(1, data.rawSignal, reshape(singlePrfData.rawSignal, 1, size(singlePrfData.rawSignal, 1), size(singlePrfData.rawSignal, 2)));
    data.isRain = cat(2, data.isRain, singlePrfData.isRain);
    data.temperature = cat(2, data.temperature, singlePrfData.temperature);
    data.azimuthAng = cat(2, data.azimuthAng, singlePrfData.azimuthAng);
    data.zenithAng = cat(2, data.zenithAng, singlePrfData.zenithAng);
    data.latitude = cat(2, data.latitude, singlePrfData.latitude);
    data.longitude = cat(2, data.longitude, singlePrfData.longitude);
    data.asl = cat(2, data.asl, singlePrfData.asl);
    data.speed = cat(2, data.speed, singlePrfData.speed);
end

end