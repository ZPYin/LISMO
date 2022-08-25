function [oz] = OverlapFunc(height, varargin)
% OVERLAPFUNC Calculate overlap function.
%
% USAGE:
%    [oz] = OverlapFunc(height)
%
% INPUTS:
%    height: numeric
%        height bins. (m)
%
% KEYWORDS:
%    channel: char
%        channel tag.
%        'far_range' or 'near_range'
%
% OUTPUTS:
%    oz: numeric
%        overlap function at each height bin.
%
% HISTORY:
%    2022-08-25: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'height', @isnumeric);
addParameter(p, 'channel', '', @ischar);

parse(p, height, varargin{:});

switch lower(p.Results.channel)
case 'far_range'
    oz = ones(size(height));
case 'near_range'
    oz = ones(size(height));
otherwise
    error('Unknown lidar channel.');
end

end