function [oData] = readVisLidarL0(file)
% READVISLIDARL0 read visibility lidar data.
%
% USAGE:
%    [oData] = readVisLidarL0(file)
%
% INPUTS:
%    file: char
%        absolute path of visibility lidar data file.
%
% OUTPUTS:
%    oData: struct
%        startTime: datenum
%        laserTemp: double
%        hRes: double
%        height: numeric (m)
%        sig: numeric ()
%
% HISTORY:
%    2024-12-17: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

oData = struct();

fid = fopen(file, 'r');

thisLine = fgetl(fid);
subStrs = strsplit(thisLine, '\t');
startTime = datenum(subStrs{3}, 'yyyy-mm-dd HH:MM:SS');
fgetl(fid);
laserTemp = str2double(fgetl(fid));
fgetl(fid);
fgetl(fid);
fgetl(fid);
thisLine = fgetl(fid);
subStrs = strsplit(thisLine, '\t');
hRes = str2double(subStrs{5});
fgetl(fid);
fgetl(fid);

data = textscan(fid, '%f%f', 'delimiter', ';');

height = data{1};
sig = data{2};

fclose(fid);

oData.startTime = startTime;
oData.laserTemp = laserTemp;
oData.height = height;
oData.sig = sig;
oData.hRes = hRes;

end