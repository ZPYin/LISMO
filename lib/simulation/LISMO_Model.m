clc;
close all;
projectDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));

%% Parameter Initialization
distArr = 7.5:15:20000;   % distance array. (m)
eleAngle = 3;   % elevation angle. (degree)
laserCW = 1030;   % central wavelength of emitting laser. (nm)
repRate = 10000;   % repitition rate.
pulseEn = 100e-3;   % pulse energy. (mJ)
dTel = 0.134;   % diameter of receiving telescope. (m)
FOV_FR = 0.3;   % field of view of far-range channel. (mrad)
FOV_NR = 0.4;   % field of near-range channel. (mrad)
etaFR = 0.08;   % quantum efficiency
etaNR = 0.08;   % quantum efficiency
PB = 0.0;   % sky spectral radiance. (W m-2 sr-1 nm-1)
darkCount = 250;   % dark counts. (s-1)
accShots = 500000;   % number of accumulation shots.
acqTime = 100;   % acquisition time. (ns)
FWHM = 3;   % FWHM of narrowband interference filter. (nm)
NIFCW = 1030;   % central wavelength of narrow band interference filter. (nm)
ND_FR = 2;   % neutral density filter. (OD)
ND_NR = 2;   % neutral density filter. (OD)
seaFogType = 'sea-fog-heavy';
saveFile = 'sea-fog-heavy-nighttime-sim.mat';
figFile = 'sea-fog-heavy-nighttime-sim.png';

%% Atmosphere Module
height = distArr * sin(eleAngle / 180 * pi);
[mBsc, mExt] = MolModel(height, laserCW, 'meteor', 'standard_atmosphere');
[aBsc, aExt] = AerModel(height, 'scene', 'marine-moderate');
[fBsc, fExt] = SeaFogModel(distArr, laserCW, 'scene', seaFogType);
tBsc = mBsc + aBsc + fBsc;
tExt = mExt + aExt + fExt;

%% Overlap
OzFR = OverlapFunc(height, 'channel', 'far_range');
OzNR = OverlapFunc(height, 'channel', 'near_range');

%% Lidar Signal Composition
% far-range
P_FR = OzFR .* 10.^(-ND_FR) * pulseEn * 1e-3 * 3e8 * pi*(dTel / 2)^2 .* tBsc .* exp(-2 * nancumsum(tExt .* [distArr(1), diff(distArr)])) ./ distArr.^2 / 2;
Ns_FR = accShots * P_FR * acqTime * 1e-9 * etaFR * laserCW * 1e-9 / (3e8 * 6.626070040e-34);
Nb_FR = accShots * 10.^(-ND_FR) * PB * pi * (FOV_FR*1e-3 / 2) ^2 * FWHM * pi * (dTel / 2)^2 * acqTime * 1e-9 * etaFR * laserCW * 1e-9 / (3e8 * 6.626070040e-34);
Nd_FR = accShots * darkCount * acqTime * 1e-9;
N_FR = Ns_FR + Nb_FR + Nd_FR;

% near-range
P_NR = OzNR .* 10.^(-ND_NR) * pulseEn * 1e-3 * 3e8 * pi*(dTel / 2)^2 .* tBsc .* exp(-2 * nancumsum(tExt .* [distArr(1), diff(distArr)])) ./ distArr.^2 / 2;
Ns_NR = accShots * P_NR * acqTime * 1e-9 * etaNR * laserCW * 1e-9 / (3e8 * 6.626070040e-34);
Nb_NR = accShots * 10.^(-ND_NR) * PB * pi * (FOV_NR*1e-3 / 2) ^2 * FWHM * pi * (dTel / 2)^2 * acqTime * 1e-9 * etaNR * laserCW * 1e-9 / (3e8 * 6.626070040e-34);
Nd_NR = accShots * darkCount * acqTime * 1e-9;
N_NR = Ns_NR + Nb_NR + Nd_NR;

%% Shot noise generation
N_FR_Poiss = sigGenWithNoise(N_FR, 0, 1, 'poisson');
N_NR_Poiss = sigGenWithNoise(N_NR, 0, 1, 'poisson');

%% data visualization

% signal profiles
subfig = subfigPos([0.08, 0.12, 0.88, 0.82], 1, 3, 0.03, 0);
figure('position', [0, 10, 750, 500], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', subfig(1, :), 'Units', 'Normalized');
sigNR = N_NR / (accShots * acqTime * 1e-3);
sigFR = N_FR / (accShots * acqTime * 1e-3);
sigFR(sigFR <= 0) = NaN;
sigNR(sigNR <= 0) = NaN;
p1 = semilogx(sigFR, distArr / 1e3, 'Color', [255, 100, 72]/255, 'LineWidth', 2, 'DisplayName', 'Far-range'); hold on;
p2 = semilogx(sigNR, distArr / 1e3, 'Color', [154, 205, 50]/255, 'LineWidth', 2, 'DisplayName', 'Near-range');

xlabel('Photon Count Rate (MHz)');
ylabel('Distance (km)');

ylim([0, 20]);
xlim([0.01, 1e6]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'XTick', logspace(-2, 6, 5), 'Box', 'on', 'FontSize', 11);

subplot('Position', subfig(2, :), 'Units', 'Normalized');
sigFR = N_FR_Poiss - Nb_FR - Nd_FR;
sigNR = N_NR_Poiss - Nb_NR - Nd_NR;
sigFR(sigFR <= 0) = NaN;
sigNR(sigNR <= 0) = NaN;
p1 = semilogx(sigFR / (accShots * acqTime * 1e-3), distArr / 1e3, 'Color', [255, 100, 72]/255, 'LineWidth', 2, 'DisplayName', 'Far-range'); hold on;
p2 = semilogx(sigNR / (accShots * acqTime * 1e-3), distArr / 1e3, 'Color', [154, 205, 50]/255, 'LineWidth', 2, 'DisplayName', 'Near-range');
p3 = semilogx((Nb_NR + Nd_NR) * ones(size(distArr)) / (accShots * acqTime * 1e-3), distArr / 1e3, 'Color', [154, 205, 50]/255, 'LineStyle', '--', 'LineWidth', 1, 'DisplayName', 'Near-range bg');
p4 = semilogx((Nb_FR + Nd_FR) * ones(size(distArr)) / (accShots * acqTime * 1e-3), distArr / 1e3, 'Color', [255, 100, 72]/255, 'LineStyle', '--', 'LineWidth', 1, 'DisplayName', 'Far-range bg');

xlabel('Signal /w noise (MHz)');
ylabel('');

ylim([0, 20]);
xlim([1e-4, 1e6]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTick', '', 'XTick', logspace(-4, 6, 6), 'Box', 'on', 'FontSize', 11);

subplot('Position', subfig(3, :), 'Units', 'Normalized');
SNRFR = (N_FR_Poiss - Nb_FR - Nd_FR) ./ sqrt(N_FR_Poiss);
SNRNR = (N_NR_Poiss - Nb_NR - Nd_NR) ./ sqrt(N_NR_Poiss);
SNRFR(SNRFR <= 0) = NaN;
SNRNR(SNRNR <= 0) = NaN;
p1 = semilogx(SNRFR, distArr / 1e3, 'Color', [255, 100, 72]/255, 'LineWidth', 2, 'DisplayName', 'Far-range'); hold on;
p2 = semilogx(SNRNR, distArr / 1e3, 'Color', [154, 205, 50]/255, 'LineWidth', 2, 'DisplayName', 'Near-range');
p5 = semilogx(3 * ones(size(distArr)), distArr / 1e3, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1, 'DisplayName', 'Detection Thresh.');

xlabel('SNR');
ylabel('');
legend([p1, p2, p3, p4, p5], 'Location', 'NorthEast');

ylim([0, 20]);
xlim([0.1, 1e6]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTick', '', 'XTick', logspace(0, 6, 4), 'Box', 'on', 'FontSize', 11);

export_fig(gcf, fullfile(projectDir, 'image', figFile), '-r300');

%% Save results
save(fullfile(projectDir, 'results', saveFile), 'eleAngle', 'distArr', 'aBsc', 'aExt', 'fBsc', 'fExt', 'mBsc', 'mExt', 'laserCW', 'repRate', 'pulseEn', 'dTel', 'FOV', 'etaFR', 'etaNR', 'PB', 'darkCount', 'accShots', 'acqTime', 'FWHM', 'NIFCW', 'ND_FR', 'ND_NR', 'seaFogType', 'N_FR_Poiss', 'N_NR_Poiss', 'OzFR', 'OzNR');