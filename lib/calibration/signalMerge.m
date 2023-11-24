function [sigGlue] = sigMerge(sigH, sigL, height, mergeRange, slope, offset)
% SIGMERGEREAL signal merge for REAL lidar data.
%
% USAGE:
%    [sigGlue] = sigMerge(sigH, sigL, height, mergeRange, slope, offset)
%
% INPUTS:
%    sigH: numeric (height x time)
%    sigL: numeric
%    height: numeric
%    mergeRange: matrix (2 x nChs)
%        signal merge range. (m)
%    slope: numeric
%        signal merge slope.
%    offset: numeric
%        signal merge offset.
%
% OUTPUTS:
%    sigGlue: numeric
%
% HISTORY:
%    2021-09-20: first edition by Zhenping
% .. Authors: - zhenping@tropos.de

height = reshape(height, size(sigH, 1), 1);
weightH = (height >= mergeRange(2)) + ...
    ((height >= mergeRange(1)) & (height < mergeRange(2))) .* ...
    ((height - mergeRange(1)) ./ (mergeRange(2) - mergeRange(1)));
weightL = 1 - weightH;

sigGlue = repmat(weightL, 1, size(sigH, 2)) .* (sigL * slope + offset) + ...
    repmat(weightH, 1, size(sigH, 2)) .* sigH;

end