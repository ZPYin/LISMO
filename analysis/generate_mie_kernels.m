clc;
close all;
projectDir = fileparts(fileparts(mfilename('fullpath')));

%% Parameter definition
wavelength = 0.21:0.01:4;
r = 0.01:0.01:100;

%% Mie scattering
qe = zeros(length(wavelength), length(r));
qs = zeros(length(wavelength), length(r));
qb = zeros(length(wavelength), length(r));
g = zeros(length(wavelength), length(r));
for iWL = 1:length(wavelength)
    fprintf('Finished %6.2f%%\n', (iWL - 1) / length(wavelength) * 100);
    [refIdxReal, refIdxImg] = refIdxWater(wavelength(iWL) * 1e3);

    for iR = 1:length(r)
        m = refIdxReal + 1i * refIdxImg;
        k0 = 2 * pi / (wavelength(iWL) * 1e3);
        a = r(iR) * 1000;

        res = Mie(m, k0 * a);

        qe(iWL, iR) = res(1);
        qs(iWL, iR) = res(2);
        qb(iWL, iR) = res(4);
        g(iWL, iR) = res(5);
    end
end

save(fullfile(projectDir, 'data', 'mie_kernel_water_droplets.mat'), 'wavelength', 'r', 'qe', 'qs', 'qb', 'g');