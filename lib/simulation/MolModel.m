function [mBsc, mExt] = MolModel(height, wavelength, varargin)
% MOLMODEL Generate molecular scattering matrix.
%
% USAGE:
%    [mBsc, mExt] = MolModel(height, wavelength)
%
% INPUTS:
%    height: numeric
%        array of altitude bins. (m)
%    wavelength: numeric
%        wavelength. (nm)
%
% KEYWORDS:
%    meteor: char
%        meteorological dataset. ('standard_atmosphere', 'radiosonde')
%
% OUTPUTS:
%    mBsc: numeric
%        molecular backscatter. (m^-1sr^-1)
%    mExt: numeric
%        molecular extinction. (m^-1sr^-1)
%
% HISTORY:
%    2022-08-19: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'height', @isnumeric);
addRequired(p, 'wavelength', @isnumeric);
addParameter(p, 'meteor', '', @ischar);

parse(p, height, wavelength, varargin{:});

%% Read meteorological data
switch lower(p.Results.meteor)
case 'standard_atmosphere'
    [alt, ~, ~, temp, pres] = atmo(max(height / 1e3) + 1, 0.03, 1);
case 'radiosonde'
otherwise
    error('Unknown meteorological data type.');
end

%% Interpolate data
presInterp = interp1(alt * 1e3, pres, height);
tempInterp = interp1(alt * 1e3, temp, height);

%% Calculate Rayleigh scattering
[mBsc, mExt] = rayleigh_scattering(wavelength, presInterp / 1e2, tempInterp, 360, 80);

end