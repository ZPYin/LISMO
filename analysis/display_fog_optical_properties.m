clc;
close all;

%% Parameter definition
wavelength = 4;
r = 0.01:0.1:60;

%% Mie scattering
qe = zeros(1, length(r));
qs = zeros(1, length(r));
[refIdxReal, refIdxImg] = refIdxWater(wavelength * 1e3);
for iR = 1:length(r)
    m = refIdxReal + 1i * refIdxImg;
    k0 = 2 * pi / (wavelength * 1e3);
    a = r(iR) * 1000;

    res = Mie(m, k0 * a);

    qe(iR) = res(1);
    qs(iR) = res(2);
end

%% Fog size dist.
fogSD = fogSD_SS('type', 'advection', 'intensity', 'moderate');

backscatter = sum(fogSD(r) .* qs .* pi .* r.^2 * 1e-9 * (r(2) - r(1)));
extinction = sum(fogSD(r) .* qe .* pi .* r.^2 * 1e-9 * (r(2) - r(1)));