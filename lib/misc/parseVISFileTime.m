function [fileTime] = parseVISFileTime(filename)
% PARSEVISFILETIME parse the file time from the filename.
%
% USAGE:
%    [fileTime] = parseVISFileTime(filename)
%
% INPUTS:
%    filename: char
%
% OUTPUTS:
%    fileTime: datenum
%
% HISTORY:
%    2024-01-10: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

filename = basename(filename);
fileTime = datenum(filename((end - 17):(end - 4)), 'yyyymmddHHMMSS');

end