function [refIdxReal, refIdxImg] = refIdxWater(inWL, varargin)
% REFIDXWATER Calculate refractive index of water at a given wavelength.
%
% USAGE:
%    [refIdxReal, refIdxImg] = refIdxWater(inWL)
%
% INPUTS:
%    inWL: numeric
%        input wavelength. (nm)
%
% KEYWORDS:
%    method: char
%        interpolation method. (default, 'linear')
%
% OUTPUTS:
%    refIdxReal: numeric
%        real part of refractive index of water.
%    refIdxImg: numeric
%        imaginary part of refractive index of water.
%
% EXAMPLE:
%
% HISTORY:
%    2022-06-29: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'inWL', @isnumeric);
addParameter(p, 'method', 'linear', @ischar);

parse(p, inWL, varargin{:});

inWL = inWL / 1e3;   % (micron)
res = load('refractive_index_water.mat');

if (inWL > res.wavelength(end)) || (inWL < res.wavelength(1))
    error('Wavelength is out of range.');
end

refIdxReal = interp1(res.wavelength, res.refractive_index_real, inWL, p.Results.method);
refIdxImg = interp1(res.wavelength, res.refractive_index_imaginary, inWL, p.Results.method);

end