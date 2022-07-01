function [fun] = fogSD_SS(varargin)
% FOGSD_SS define fog size distribution based on Silverman-Sprague model.
%
% USAGE:
%    [fun] = fogSD_SS()
%
% INPUTS:
%
% KEYWORDS:
%    type: char
%        'advection' or 'radiation' fog.
%    intensity: char
%        'moderate', 'heavy'
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

addParameter(p, 'type', 'advection', @ischar);
addParameter(p, 'intensity', 'moderate', @ischar);

parse(p, varargin{:});

switch lower([p.Results.type, '+', p.Results.intensity])
case 'advection+heavy'
    a = 3;
    b = 1;
    c = 0.027e6;
    d = 0.3;
case 'advection+moderate'
    a = 3;
    b = 1;
    c = 0.06592e6;
    d = 0.375;
case 'radiation+heavy'
    a = 6;
    b = 1;
    c = 2.37305e6;
    d = 1.5;
case 'radiation+moderate'
    a = 6;
    b = 1;
    c = 607.5e6;
    d = 3.0;
otherwise
    error('Unknown fog type.');
end

fun = @(r) c .* r.^a .* exp(-d .* r.^b);

end