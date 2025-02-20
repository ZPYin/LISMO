function [oData] = readVIS_ZKGB(dataFile, varargin)
% READVIS_ZKGB read visibility lidar data from ZKGB.
%
% USAGE:
%    [oData] = readVIS_ZKGB(dataFile)
%
% INPUTS:
%    dataFile: char
%        The path of the data file.
%
% KEYWORDS:
%    debug: logical
%
% OUTPUTS:
%    oData: struct
%        nBins: int
%        range: (m)
%        angle: (degree)
%        signal: angle x range
%
% HISTORY:
%    2025-02-08: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'dataFile', @ischar);
addParameter(p, 'debug', false, @islogical);

parse(p, dataFile, varargin{:});

if p.Results.debug
    fprintf('Reading data from %s\n', dataFile);
end

thisDataFilename = basename(dataFile);
startTime = datenum(thisDataFilename(1:14), 'yyyymmddHHMMSS');

oData = struct();
thisMat = readmatrix(dataFile, 'range', [1 1]);
oData.nBins = size(thisMat, 2) - 1;
oData.range = thisMat(1, 2:end);
oData.angle = thisMat(2:end, 1);
oData.signal = thisMat(2:end, 2:end);
oData.mTime = startTime + datenum(0, 1, 0, 0, 20 / 90 * (0:(size(oData.signal, 1) - 1)), 0);

end