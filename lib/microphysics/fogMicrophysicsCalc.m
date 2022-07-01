function [res] = fogMicrophysicsCalc(fun,varargin)
% FOGMICROPHYSICSCALC Calculate microphysics of fog.
%
% USAGE:
%    [res] = fogMicrophysicsCalc(fun)
%
% INPUTS:
%    fun: function handle
%        function handle of fog size distribution. (m-3*um-1)
%
% KEYWORDS:
%    xmin: numeric
%        minimum integral limit. (um)
%    xmax: numeric
%        maximum integral limit. (um)
%
% OUTPUTS:
%    res: struct
%        number_concentration
%        mean_radius
%        effective_radius
%
% EXAMPLE:
%
% HISTORY:
%    2022-06-30: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'fun', @(x)isa(x, 'function_handle'));
addParameter(p, 'xmin', 0.01, @isnumeric);
addParameter(p, 'xmax', 60, @isnumeric);

parse(p, fun, varargin{:});

res = struct();
res.number_concentration = integral(fun, p.Results.xmin, p.Results.xmax);
res.mean_radius = integral(@(r)fun(r).*r, p.Results.xmin, p.Results.xmax) / res.number_concentration;
res.effective_radius = integral(@(r)fun(r) .* r.^3, p.Results.xmin, p.Results.xmax) / integral(@(r)fun(r) .* r.^2, p.Results.xmin, p.Results.xmax);

end