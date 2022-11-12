function [visibility, visibilityErr] = ext2vis(tExt, varargin)
% EXT2VIS convert extinction coefficient to visibility.
%
% USAGE:
%    [visibility] = ext2vis(tExt)
%
% INPUTS:
%    tExt: array
%        total extinction coefficient. (m^-1)
%
% KEYWORDS:
%    method: char
%        visibility calculation method. (default: 'wmo')
%    ExtinctionError: array
%        error of extinction coefficient. (m^-1)
%
% OUTPUTS:
%    visibility: array
%        visibility. (m)
%    visibilityErr: array
%        error of visibility. (m)
%
% HISTORY:
%    2022-11-09: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'tExt', @isnumeric);
addParameter(p, 'method', 'wmo', @ischar);
addParameter(p, 'ExtinctionError', [], @isnumeric);

parse(p, tExt, varargin{:});

switch lower(p.Results.method)
case 'wmo'
    visibility = 3 ./ tExt;

    if ~ isempty(p.Results.ExtinctionError)
        visibilityErr = sqrt(1 ./ tExt.^2 * p.Results.ExtinctionError);
    end

otherwise
    error('Unknown visibility calculation method.');
end

end