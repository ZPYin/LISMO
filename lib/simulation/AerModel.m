function [aBsc, aExt] = AerModel(height, varargin)
% AERMODEL Generate aerosol optical properties.
%
% USAGE:
%    [aBsc, aExt] = AerModel(height)
%
% INPUTS:
%    height: numeric
%        array of altitude bins. (m)
%
% KEYWORDS:
%    scene: char
%        aerosol conditions.
%        'marine-weak' (default)
%        'marine-moderate'
%        'marine-heavy'
%
% OUTPUTS:
%    aBsc: numeric
%        aerosol backscatter. (m^-1sr^-1)
%    aExt: numeric
%        aerosol extinction. (m^-1sr^-1)
%
% HISTORY:
%    2022-08-19: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'height', @isnumeric);
addParameter(p, 'scene', 'marine-weak', @ischar);

parse(p, height, varargin{:});

%% Read meteorological data
switch lower(p.Results.scene)
case 'marine-weak'
    layerBase = 0;   % distance of sea fog base. (m)
    layerTop = 1500;   % distance of sea fog top. (m)
    lr = 20;   % lidar ratio of sea fog. (sr)
    aExtMax = 50e-6;   % maximum extinction at 532 nm. (m-1)
    aExtFunc = @(h) (-aExtMax / (layerTop - layerBase).^2 .* (h - layerBase).^2 + aExtMax) .* ((h >= layerBase) & (h <= layerTop));
    aBscFunc = @(h) aExtFunc(h) / lr;
case 'marine-moderate'
    layerBase = 0;   % distance of sea fog base. (m)
    layerTop = 1500;   % distance of sea fog top. (m)
    lr = 20;   % lidar ratio of sea fog. (sr)
    aExtMax = 100e-6;   % maximum extinction at 532 nm. (m-1)
    aExtFunc = @(h) (-aExtMax / (layerTop - layerBase).^2 .* (h - layerBase).^2 + aExtMax) .* ((h >= layerBase) & (h <= layerTop));
    aBscFunc = @(h) aExtFunc(h) / lr;
case 'marine-heavy'
    layerBase = 0;   % distance of sea fog base. (m)
    layerTop = 1500;   % distance of sea fog top. (m)
    lr = 20;   % lidar ratio of sea fog. (sr)
    aExtMax = 200e-6;   % maximum extinction at 532 nm. (m-1)
    aExtFunc = @(h) (-aExtMax / (layerTop - layerBase).^2 .* (h - layerBase).^2 + aExtMax) .* ((h >= layerBase) & (h <= layerTop));
    aBscFunc = @(h) aExtFunc(h) / lr;
otherwise
    error('Unknown meteorological data type.');
end

aBsc = aBscFunc(height);
aExt = aExtFunc(height);

end