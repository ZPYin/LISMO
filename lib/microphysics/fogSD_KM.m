function [fun] = fogSD_KM(wc, varargin)
% FOGSD_KM define fog size distribution based on Khragian-Mazin model.
%
% USAGE:
%    [fun] = fogSD_KM(wc)
%
% INPUTS:
%    wc: numeric
%        liquid water content. (g*m^-3)
%
% KEYWORDS:
%    type: char
%        'advection' or 'radiation' fog.
%
% OUTPUTS:
%    fun: handle
%        size distribution function handle. (m^-3*um^-1)
%
% HISTORY:
%    2022-06-30: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'wc', @isnumeric);
addParameter(p, 'type', 'advection', @ischar);

parse(p, wc, varargin{:});

switch lower(p.Results.type)
case 'advection'
    fun = @(r) 3.73 * 1e5 * wc.^-0.804 * r.^2 .* exp(-0.2392 * wc.^-0.301 .* r);
case 'radiation'
    fun = @(r) 5.4 * 1e7 * wc.^-1.104 * r.^2 .* exp(-0.5477 * wc.^-0.351 .* r);
otherwise
    error('Unknown fog type.');
end

end