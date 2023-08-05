function [data] = readV35(file, varargin)
% READV35 read visibility data.
%
% USAGE:
%    [data] = readV35(file)
%
% INPUTS:
%    file: char
%        absolute path of data file.
%
% OUTPUTS:
%    data: struct
%        mTime: numeric
%            measurement time in datanum.
%        vis1min: numeric
%            visibility of 1 min average.
%        vis10min: numeric
%            visibility of 10 min average.
%        temp: numeric
%            ambient temperature. (C)
%
% HISTORY:
%    2022-07-14: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'file', @ischar);

parse(p, file, varargin{:});

if exist(file, 'file') ~= 2
    error('file does not exist!');
end

fid = fopen(file, 'r');

dataInfo = textscan(fid, '%s%f%f%f', 'HeaderLines', 1, 'Delimiter', ',');

fclose(fid);

mTime = zeros(1, length(dataInfo{1}));
vis1min = NaN(1, length(dataInfo{1}));
vis10min = NaN(1, length(dataInfo{1}));
temp = NaN(1, length(dataInfo{1}));
for iR = 1:length(dataInfo{1})
    try
        mTime(iR) = datenum(dataInfo{1}{iR}, 'yyyy-mm-dd HH:MM:SS');
        vis1min(iR) = dataInfo{2}(iR);
        vis10min(iR) = dataInfo{3}(iR);
        temp(iR) = dataInfo{4}(iR);
    catch
        warning('failure in parsing data entry.');
    end
end

data = struct();
data.mTime = mTime;
data.vis1min = vis1min;
data.vis10min = vis10min;
data.temp = temp;

end