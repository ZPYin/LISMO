 function [tExt] = vis2ext(visibility, varargin)
% VIS2EXT convert extinction coefficient to visibility.
%
% USAGE:
%    [visibility] = vis2ext(tExt)
%
% INPUTS:
%    visibility: array
%        visibility. (m)
%
% KEYWORDS:
%    method: char
%        visibility calculation method. (default: 'wmo', 'mor')
%
% OUTPUTS:
%    tExt: array
%        total extinction coefficient. (m^-1)
%
% HISTORY:
%    2022-11-22: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'visibility', @isnumeric);
addParameter(p, 'method', 'wmo', @ischar);

parse(p, visibility, varargin{:});

switch lower(p.Results.method)
case 'wmo'
    tExt = -log(0.02) ./ visibility;

case 'mor'
    tExt = 3 ./ visibility;

otherwise
    error('Unknown visibility calculation method.');
end

end