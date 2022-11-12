function [sfBsc, sfExt] = SeaFogModel(distance, wavelength, varargin)
% SEAFOGMODEL Generate sea fog optical properties.
%
% USAGE:
%    [sfBsc, sfExt] = SeaFogModel(distance, wavelength)
%
% INPUTS:
%    distance: numeric
%        array of range bins. (m)
%
% KEYWORDS:
%    scene: char
%        sea fog conditions.
%        'sea-fog-weak'
%        'sea-fog-moderate'
%        'sea-fog-heavy'
%        'sea-fog-none'
%    seaFogType: char
%        sea fog type.
%        'moderate-radiation-fog' (default)
%        'heavy-radiation-fog'
%        'moderate-advection-fog'
%        'heavy-advection-fog'
%    distLayerFront: numeric
%        distance of front edge of sea fog to the lidar. (m)
%    distLayerBack: numeric
%        distance of back edge of sea fog to the lidar. (m)
%
% OUTPUTS:
%    sfBsc: numeric
%        sea fog backscatter. (m^-1sr^-1)
%    sfExt: numeric
%        sea fog extinction. (m^-1sr^-1)
%
% HISTORY:
%    2022-08-19: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'distance', @isnumeric);
addRequired(p, 'wavelength', @isnumeric);
addParameter(p, 'scene', 'sea-fog-moderate', @ischar);
addParameter(p, 'seaFogType', 'moderate-radiation-fog', @ischar);
addParameter(p, 'distLayerFront', 4000, @isnumeric);
addParameter(p, 'distLayerBack', 5000, @isnumeric);

parse(p, distance, wavelength, varargin{:});

%% Load sea fog optical properties
projectDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
sfOpt = load(fullfile(projectDir, 'data', 'sea-fog-lidar-ratio.mat'));
switch lower(p.Results.seaFogType)
case 'moderate-radiation-fog'
    lr = interp1(sfOpt.wavelength * 1e3, sfOpt.lrRadModerate, wavelength);
case 'heavy-radiation-fog'
    lr = interp1(sfOpt.wavelength * 1e3, sfOpt.lrRadHeavy, wavelength);
case 'moderate-advection-fog'
    lr = interp1(sfOpt.wavelength * 1e3, sfOpt.lrAdvModerate, wavelength);
case 'heavy-advection-fog'
    lr = interp1(sfOpt.wavelength * 1e3, sfOpt.lrAdvHeavy, wavelength);
otherwise
    error('Unknown sea fog type.');
end

%% Read meteorological data
switch lower(p.Results.scene)
case 'sea-fog-weak'
    sfExtMax = 500e-6;   % maximum extinction. (m-1)
    sfExtFunc = @(r) (sfExtMax / (p.Results.distLayerBack - p.Results.distLayerFront)^2 * 4 .* (r - p.Results.distLayerFront) .* (p.Results.distLayerBack - r)) .* ((r >= p.Results.distLayerFront) & (r <= p.Results.distLayerBack));
    sfBscFunc = @(r) sfExtFunc(r) / lr;
case 'sea-fog-moderate'
    sfExtMax = 1000e-6;   % maximum extinction. (m-1)
    sfExtFunc = @(r) (sfExtMax / (p.Results.distLayerBack - p.Results.distLayerFront)^2 * 4 .* (r - p.Results.distLayerFront) .* (p.Results.distLayerBack - r)) .* ((r >= p.Results.distLayerFront) & (r <= p.Results.distLayerBack));
    sfBscFunc = @(r) sfExtFunc(r) / lr;
case 'sea-fog-heavy'
    sfExtMax = 2000e-6;   % maximum extinction. (m-1)
    sfExtFunc = @(r) (sfExtMax / (p.Results.distLayerBack - p.Results.distLayerFront)^2 * 4 .* (r - p.Results.distLayerFront) .* (p.Results.distLayerBack - r)) .* ((r >= p.Results.distLayerFront) & (r <= p.Results.distLayerBack));
    sfBscFunc = @(r) sfExtFunc(r) / lr;
case 'sea-fog-none'
    sfExtFunc = @(r) zeros(size(r));
    sfBscFunc = @(r) zeros(size(r));
otherwise
    error('Unknown meteorological data type.');
end

sfBsc = sfBscFunc(distance);
sfExt = sfExtFunc(distance);

end