function [data] = readSingleVIS(filename, varargin)
% READSINGLEVIS Read fog lidar from WXZK.
%
% USAGE:
%    [data] = readSingleVIS(filename)
%
% INPUTS:
%    filename: char
%        absolute path of data file.
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

addRequired(p, 'filename', @ischar);
addParameter(p, 'debug', false, @islogical);

parse(p, filename, varargin{:});

data = struct();
nChs = 4;

if exist(filename, 'file') ~= 2
    warning('file does not exist %s', filename);
    return;
end

try

    fid = fopen(filename, 'r');

    % first line is the actual filename of the profile.
    % i.e.: WX20220818135403.VIS
    fgetl(fid);

    % second line contains system information for the current profile.
    % i.e.:  WX 08/18/2022 13:54:03 08/18/2022 13:54:03 E:119.369485 N:25.3922722 0 0 0 0 0 3 0 0 
    line2 = fgetl(fid);
    subStrs = strsplit(strip(line2), ' ');
    mTimeStart = datenum([subStrs{2}, subStrs{3}], 'mm/dd/yyyyHH:MM:SS');
    mTimeStop = datenum([subStrs{4}, subStrs{5}], 'mm/dd/yyyyHH:MM:SS');
    longitude = str2double(subStrs{6}(3:end));
    latitude = str2double(subStrs{7}(3:end));
    asl = str2double(subStrs{8});
    speed = str2double(subStrs{9});
    azimuthAng = str2double(subStrs{12});
    zenithAng = 90.0 - str2double(subStrs{13});
    temperature = str2double(subStrs{14});
    isRain = strcmp(subStrs{15}, '1');

    % third line contains laser status information.
    % i.e.: 5000 5000 
    fgetl(fid);

    % fourth line (to the end of header) contains channel information.
    % 0 3000 7.5 1064.n 64 60000 1000 
    % 1 3000 7.5 1064.p 64 60000 1000 
    % 2 3000 7.5 1064.s 64 60000 1000 
    % 3 3000 7.5 1064.min 64 60000 1000 
    line4 = fgetl(fid);
    subStrs = strsplit(strip(line4), ' ');
    nBins = str2double(subStrs{2});
    hRes = str2double(subStrs{3});
    nShots = str2double(subStrs{6});

    for iCh = 1:(nChs - 1)
        fgetl(fid);
    end

    fread(fid, 2, 'uint8');

    %% read data matrix
    chSignal = NaN(nChs, nBins);
    for iCh = 1:nChs
        chSignal(iCh, :) = fread(fid, nBins, 'uint64');
    end

    fclose(fid);

catch ME
    if p.Results.debug
        fprintf('Cracked file: %s\n', filename);
    end

    fclose(fid);
    rethrow(ME)
end

data.hRes = hRes;
data.nShots = nShots;
data.nBins = nBins;
data.startTime = mTimeStart;
data.stopTime = mTimeStop;
data.rawSignal = chSignal;
data.isRain = isRain;
data.temperature = temperature;
data.azimuthAng = azimuthAng;
data.zenithAng = zenithAng;
data.latitude = latitude;
data.longitude = longitude;
data.asl = asl;
data.speed = speed;

end