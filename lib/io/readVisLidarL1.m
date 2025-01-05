function [data] = readVisLidarL1(file)
% READVISLIDARL1 read product from visiblity lidar level 1 file.
%
% USAGE:
%    [data] = readVisLidarL1(file)
%
% INPUTS:
%    file: char
%        absolute path of data file.
%
% OUTPUTS:
%    data: struct
%        height: range bin. (m)
%        extinction: extinction. (km-1)
%        vis: visibility. (km)
%        pm: particulate matter. (ug/m3)
%        mTime: measurement time.
%
% HISTORY:
%    2024-12-19: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

fid = fopen(file, 'r');

thisLine = fgetl(fid);
subStrs = strsplit(thisLine, '\t');

tmp = textscan(fid, '%f%f%f%f', 'Delimiter', ';', 'HeaderLines', 14);

fclose(fid);

data.height = tmp{1};
data.extinction = tmp{2};
data.vis = tmp{3};
data.pm = tmp{4};
data.mTime = datenum(subStrs{2}, 'yyyy-mm-dd HH:MM:SS');

end