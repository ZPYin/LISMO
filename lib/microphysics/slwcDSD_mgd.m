function [fun] = slwcDSD_mgd(varargin)
% SLWCDSD_MGD define stratiform low-level cloud droplet size distribution based on modified gamma distribution.
%
% USAGE:
%    [fun] = slwcDSD_mgd()
%
% KEYWORDS:
%    type: char
%        'marine' or 'continental' clouds.
%    N0: numeric
%        total number concentration (m^-3). Required when type is 'input'.
%    gamma0: numeric
%        shape parameter (dimensionless). Required when type is 'input'.
%    Rm: numeric
%        mode radius (um). Required when type is 'input'.
%
% OUTPUTS:
%    fun: handle
%        size distribution function handle. (m^-3*um^-1)
%
% References:
%   Miles NL, Verlinde J, Clothiaux EE. Cloud droplet size distributions in low-level stratiform clouds. Journal of the atmospheric sciences. 2000;57(2):295-311.
%
% HISTORY:
%    2025-08-09: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn


p = inputParser;
p.KeepUnmatched = true;

addParameter(p, 'type', 'input', @ischar);
addParameter(p, 'N0', [], @isnumeric);
addParameter(p, 'gamma0', [], @isnumeric);
addParameter(p, 'Rm', [], @isnumeric);

parse(p, varargin{:});

switch lower(p.Results.type)
case 'marine'
    N0 = 74 * 1e6;   % total number concentration (m^-3)
    gamma0 = 8.6;   % shape parameter (dimensionless)
    Rm = 2.7 / 2;   % mode radius (um); Reff / (gamma + 2)
case 'continental'
    N0 = 288 * 1e6;
    gamma0 = 8.7;
    Rm = 1.3 / 2;
case 'input'
    N0 = p.Results.N0;
    gamma0 = p.Results.gamma0;
    Rm = p.Results.Rm;
otherwise
    error('Unknown cloud type.');
end

fun = @(r) N0 / Rm ./ gamma(gamma0) .* (r ./ Rm) .^ (gamma0 - 1) .* exp(- r ./ Rm);

end