function [dataSim] = LISMO_Model(distArr, varargin)
% LISMO_MODEL Simulator for sea-fog lidar backscatter signal.
%
% USAGE:
%    [dataSim] = LISMO_Model(distArr)
%
% INPUTS:
%    distArr: array
%        distance array. (m)
%
% KEYWORDS:
%    eleAngle: double
%        elevation angle. (default: 3 degree)
%    laserWL: double
%        emitting laser wavelength. (default: 1030 nm)
%    pulseEn: double
%        pulse energy. (default: 0.1 mJ)
%    dTel: double
%        diameter of receiving telescope. (default: 0.134 m)
%    FOV_FR: double
%        field of view of far-range detection channel. (default: 0.3 mrad)
%    FOV_FR: double
%        field of view of near-range detection channel. (default: 0.4 mrad)
%    etaFR: double
%        quantum efficiency of far-range detector. (default: 0.08)
%    etaNR: double
%        quantum efficiency of near-range detector. (default: 0.08)
%    PB: double
%        sky spectral radiance. (default: 0.0 W m^-2 nm^-1 sr^-1)
%    darkCountFR: double
%        dark count rate for far-range detector. (default: 250 s-1)
%    darkCountNR: double
%        dark count rate for near-range detector. (default: 250 s-1)
%    accShots: double
%        number of accumulation shots. (default: 500000)
%    acqTime: double
%        integration time for single range bin. (default: 100 ns)
%    FWHM: double
%        full width at half maximum of narrowband filter. (default: 3 nm)
%    ND_FR: double
%        optical depth of neutral density filter at far-range channel. (default: 2)
%    ND_NR: double
%        optical depth of neutral density filter at near-range channel. (default: 2)
%    tBsc: array
%        total backscatter coefficient. (m^-1sr^-1)
%    tExt: array
%        total extinction coefficient. (m^-1)
%    visible: char
%        whether to display the signal profile. ('on' (default) | 'off')
%
% OUTPUTS:
%    dataSim: struct
%        distArr: distance array. (m)
%        N_FR_Poiss: far-range signal with Poisson noise.
%        N_FR: far-range siganl without Poisson noise.
%        Nb_FR: far-range background signal.
%        Ns_FR: far-range backscattering signal.
%        N_NR_Poiss: near-range signal with Poisson noise.
%        N_NR: near-range siganl without Poisson noise.
%        Nb_NR: near-range background signal.
%        Ns_NR: near-range backscattering signal.
%
% HISTORY:
%    2022-11-09: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'distArr', @isnumeric);
addParameter(p, 'eleAngle', 3, @isnumeric);
addParameter(p, 'laserWL', 1030, @isnumeric);
addParameter(p, 'pulseEn', 100e-3, @isnumeric);
addParameter(p, 'dTel', 0.134, @isnumeric);
addParameter(p, 'FOV_FR', 0.3, @isnumeric);
addParameter(p, 'FOV_NR', 0.4, @isnumeric);
addParameter(p, 'etaFR', 0.08, @isnumeric);
addParameter(p, 'etaNR', 0.08, @isnumeric);
addParameter(p, 'PB', 0.0, @isnumeric);
addParameter(p, 'darkCountFR', 250, @isnumeric);
addParameter(p, 'darkCountNR', 250, @isnumeric);
addParameter(p, 'accShots', 500000, @isnumeric);
addParameter(p, 'acqTime', 100, @isnumeric);
addParameter(p, 'FWHM', 3, @isnumeric);
addParameter(p, 'ND_FR', 2, @isnumeric);
addParameter(p, 'ND_NR', 2, @isnumeric);
addParameter(p, 'tBsc', [], @isnumeric);
addParameter(p, 'tExt', [], @isnumeric);
addParameter(p, 'visible', 'on', @ischar);

parse(p, distArr, varargin{:});

height = distArr * sin(p.Results.eleAngle / 180 * pi);

%% Overlap
OzFR = OverlapFunc(height, 'channel', 'far_range');
OzNR = OverlapFunc(height, 'channel', 'near_range');

%% Lidar Signal Composition
% far-range
P_FR = OzFR .* 10.^(-p.Results.ND_FR) * p.Results.pulseEn * 1e-3 * 3e8 * pi*(p.Results.dTel / 2)^2 .* p.Results.tBsc .* exp(-2 * nancumsum(p.Results.tExt .* [distArr(1), diff(distArr)])) ./ distArr.^2 / 2;
Ns_FR = p.Results.accShots * P_FR * p.Results.acqTime * 1e-9 * p.Results.etaFR * p.Results.laserWL * 1e-9 / (3e8 * 6.626070040e-34);
Nb_FR = p.Results.accShots * 10.^(-p.Results.ND_FR) * p.Results.PB * pi * (p.Results.FOV_FR*1e-3 / 2) ^2 * p.Results.FWHM * pi * (p.Results.dTel / 2)^2 * p.Results.acqTime * 1e-9 * p.Results.etaFR * p.Results.laserWL * 1e-9 / (3e8 * 6.626070040e-34);
Nd_FR = p.Results.accShots * p.Results.darkCountFR * p.Results.acqTime * 1e-9;
N_FR = Ns_FR + Nb_FR + Nd_FR;

% near-range
P_NR = OzNR .* 10.^(-p.Results.ND_NR) * p.Results.pulseEn * 1e-3 * 3e8 * pi*(p.Results.dTel / 2)^2 .* p.Results.tBsc .* exp(-2 * nancumsum(p.Results.tExt .* [distArr(1), diff(distArr)])) ./ distArr.^2 / 2;
Ns_NR = p.Results.accShots * P_NR * p.Results.acqTime * 1e-9 * p.Results.etaNR * p.Results.laserWL * 1e-9 / (3e8 * 6.626070040e-34);
Nb_NR = p.Results.accShots * 10.^(-p.Results.ND_NR) * p.Results.PB * pi * (p.Results.FOV_NR*1e-3 / 2) ^2 * p.Results.FWHM * pi * (p.Results.dTel / 2)^2 * p.Results.acqTime * 1e-9 * p.Results.etaNR * p.Results.laserWL * 1e-9 / (3e8 * 6.626070040e-34);
Nd_NR = p.Results.accShots * p.Results.darkCountNR * p.Results.acqTime * 1e-9;
N_NR = Ns_NR + Nb_NR + Nd_NR;

%% Shot noise generation
N_FR_Poiss = sigGenWithNoise(N_FR, 0, 1, 'poisson');
N_NR_Poiss = sigGenWithNoise(N_NR, 0, 1, 'poisson');

dataSim = struct();
dataSim.distArr = distArr;
dataSim.N_FR_Poiss = N_FR_Poiss;
dataSim.N_FR = N_FR;
dataSim.Nb_FR = Nb_FR;
dataSim.Ns_FR = Ns_FR;
dataSim.Nd_FR = Nd_FR;
dataSim.N_NR_Poiss = N_NR_Poiss;
dataSim.N_NR = N_NR;
dataSim.Nb_NR = Nb_NR;
dataSim.Nd_NR = Nd_NR;
dataSim.Ns_NR = Ns_NR;

%% data visualization

% signal profiles
if strcmpi(p.Results.visible, 'on')
    subfig = subfigPos([0.08, 0.12, 0.88, 0.82], 1, 3, 0.03, 0);
    figure('position', [0, 10, 750, 500], 'Units', 'Pixels', 'Color', 'w');

    subplot('Position', subfig(1, :), 'Units', 'Normalized');
    sigNR = N_NR / (p.Results.accShots * p.Results.acqTime * 1e-3);
    sigFR = N_FR / (p.Results.accShots * p.Results.acqTime * 1e-3);
    sigFR(sigFR <= 0) = NaN;
    sigNR(sigNR <= 0) = NaN;
    semilogx(sigFR, distArr / 1e3, 'Color', [255, 100, 72]/255, 'LineWidth', 2, 'DisplayName', 'Far-range'); hold on;
    semilogx(sigNR, distArr / 1e3, 'Color', [154, 205, 50]/255, 'LineWidth', 2, 'DisplayName', 'Near-range');

    xlabel('Photon Count Rate (MHz)');
    ylabel('Distance (km)');

    ylim([0, 20]);
    xlim([1e-4, 1e6]);

    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'XTick', logspace(-2, 6, 5), 'Box', 'on', 'FontSize', 11);

    subplot('Position', subfig(2, :), 'Units', 'Normalized');
    sigFR = N_FR_Poiss - Nb_FR - Nd_FR;
    sigNR = N_NR_Poiss - Nb_NR - Nd_NR;
    sigFR(sigFR <= 0) = NaN;
    sigNR(sigNR <= 0) = NaN;
    semilogx(sigFR / (p.Results.accShots * p.Results.acqTime * 1e-3), distArr / 1e3, 'Color', [255, 100, 72]/255, 'LineWidth', 2, 'DisplayName', 'Far-range'); hold on;
    semilogx(sigNR / (p.Results.accShots * p.Results.acqTime * 1e-3), distArr / 1e3, 'Color', [154, 205, 50]/255, 'LineWidth', 2, 'DisplayName', 'Near-range');
    p3 = semilogx((Nb_NR + Nd_NR) * ones(size(distArr)) / (p.Results.accShots * p.Results.acqTime * 1e-3), distArr / 1e3, 'Color', [154, 205, 50]/255, 'LineStyle', '--', 'LineWidth', 1, 'DisplayName', 'Near-range bg');
    p4 = semilogx((Nb_FR + Nd_FR) * ones(size(distArr)) / (p.Results.accShots * p.Results.acqTime * 1e-3), distArr / 1e3, 'Color', [255, 100, 72]/255, 'LineStyle', '--', 'LineWidth', 1, 'DisplayName', 'Far-range bg');

    xlabel('Signal /w noise (MHz)');
    ylabel('');

    ylim([0, 20]);
    xlim([1e-4, 1e6]);

    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'XTick', logspace(-4, 6, 6), 'Box', 'on', 'FontSize', 11);

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

    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'XTick', logspace(0, 6, 4), 'Box', 'on', 'FontSize', 11);
end

end