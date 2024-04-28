function [aBsc, aExt, aBscStd, aExtStd] = extRet_Holger(r, sig, varargin)
% EXTRET_HOLGER Extinction retrieval with Quasi-Retrieval (Holger Baars, AMT, 2017)
%
% USAGE:
%    [aBsc, aExt, aBscStd, aExtStd] = extRet_Holger(r, sig)
%
% INPUTS:
%    r: numeric
%        distance of each range bin. (m)
%    sig: numeric
%        (background removed) receiving signal.
%
% KEYWORDS:
%    calibration_constant: numeric
%        lidar calibration constant. (default: 1)
%    fullOverlapR: numeric
%        minimum range distance with full overlap. (m)
%    elevation_angle: numeric
%        elevation angle. (degree, default: 0)
%    emit_wavelength: numeric
%        wavelength of emitting laser pulse. (nm)
%
% OUTPUTS:
%    aBsc: numeric
%        aerosol backscatter. (m-1sr-1)
%    aExt: numeric
%        aerosol extinction. (m-1)
%    aBscStd: numeric
%        uncertainty of aerosol backscatter. (m-1sr-1)
%    aExtStd: numeric
%        uncertainty of aerosol extinction. (m-1)
%
% HISTORY:
%    2023-04-25: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'r', @isnumeric);
addRequired(p, 'sig', @isnumeric);
addParameter(p, 'calibration_constant', 1, @isnumeric);
addParameter(p, 'fullOverlapR', 0, @isnumeric);
addParameter(p, 'elevation_angle', 0, @isnumeric);
addParameter(p, 'emit_wavelength', 1064, @isnumeric);

parse(p, r, sig, varargin{:});

%% Range Corrected signal
rcs = sig .* r.^2;

alt = r .* sin(p.Results.elevation_angle / 180 * pi);
LRaer = 50;

%% Signal calibration
attn = rcs ./ p.Results.calibration_constant;
diffR = [r(1), diff(r)];

%% Remove Blind Zone
isInBlindZone = (r <= p.Results.fullOverlapR);
attn(isInBlindZone) = NaN;

%% Obtain Rayleigh Scattering Properties
[temperature, pressure, ~, ~] = read_meteordata(0, zeros(size(r)), ...
    'meteor_data', 'standard_atmosphere', ...
    'station', 'xxx');
[mBsc, mExt] = rayleigh_scattering(p.Results.emit_wavelength, pressure, temperature + 273.14, 380, 60);

%% Extinction Profile Retrieval
mol_att = exp(- nancumsum(mExt .* diffR, 2));
quasi_par_ext = zeros(size(mBsc));

for iLoop = 1:5
    quasi_par_att = exp(-nancumsum(quasi_par_ext .* diffR, 2));
    quasi_par_bsc = attn ./ (mol_att .* quasi_par_att).^2 - mBsc;
    quasi_par_bsc(quasi_par_bsc <= 0) = 0;
    quasi_par_ext = quasi_par_bsc .* LRaer;
end

aBsc = quasi_par_att;
aExt = quasi_par_ext;
aBscStd = zeros(size(aBsc));
aExtStd = zeros(size(aExt));

end