%------------------------------------------------------------------------------%
% Simulate sea-fog lidar signal under clear-sky and sea-fog cover.
% Simulate forward/backward retrievals
% Compute errors of forward/backward retrievals
% Data visualization
%------------------------------------------------------------------------------%

global LISMO_VARS;
close all;
clc;

%% Parameter Definition
distArr = 7.5:15:20000;
eleAngle = 3;
laserWL = 1030;
hRange = [0, 10];
runSimulation = true;
savePath = fullfile(LISMO_VARS.projectDir, 'image');

%% Signal Simulation
if runSimulation

    if (~ exist(savePath, 'dir'))
        mkdir(savePath);
    end

    % clear sky
    height = distArr * sin(eleAngle / 180 * pi);
    seaFogType = 'sea-fog-none';
    [mBsc, mExt] = MolModel(height, laserWL, 'meteor', 'standard_atmosphere');
    [aBsc, aExt] = AerModel(height, 'scene', 'marine-moderate');
    [fBsc, fExt] = SeaFogModel(distArr, laserWL, 'scene', seaFogType);
    dataSim = LISMO_Model(distArr, 'tBsc', mBsc + aBsc + fBsc, 'tExt', mExt + aExt + fExt, 'eleAngle', eleAngle, 'laserWL', laserWL, 'visible', 'off');
    save(fullfile(savePath, sprintf('%s.mat', seaFogType)), 'distArr', 'eleAngle', 'aBsc', 'aExt', 'fBsc', 'fExt', 'mBsc', 'mExt', 'dataSim');

    % sea-fog-weak
    seaFogType = 'sea-fog-weak';
    [mBsc, mExt] = MolModel(height, laserWL, 'meteor', 'standard_atmosphere');
    [aBsc, aExt] = AerModel(height, 'scene', 'marine-moderate');
    [fBsc, fExt] = SeaFogModel(distArr, laserWL, 'scene', seaFogType);
    dataSim = LISMO_Model(distArr, 'tBsc', mBsc + aBsc + fBsc, 'tExt', mExt + aExt + fExt, 'eleAngle', eleAngle, 'laserWL', laserWL, 'visible', 'off');
    save(fullfile(savePath, sprintf('%s.mat', seaFogType)), 'distArr', 'eleAngle', 'aBsc', 'aExt', 'fBsc', 'fExt', 'mBsc', 'mExt', 'dataSim');

    % sea-fog-moderate
    seaFogType = 'sea-fog-moderate';
    [mBsc, mExt] = MolModel(height, laserWL, 'meteor', 'standard_atmosphere');
    [aBsc, aExt] = AerModel(height, 'scene', 'marine-moderate');
    [fBsc, fExt] = SeaFogModel(distArr, laserWL, 'scene', seaFogType);
    dataSim = LISMO_Model(distArr, 'tBsc', mBsc + aBsc + fBsc, 'tExt', mExt + aExt + fExt, 'eleAngle', eleAngle, 'laserWL', laserWL, 'visible', 'off');
    save(fullfile(savePath, sprintf('%s.mat', seaFogType)), 'distArr', 'eleAngle', 'aBsc', 'aExt', 'fBsc', 'fExt', 'mBsc', 'mExt', 'dataSim');

    % sea-fog-heavy
    seaFogType = 'sea-fog-heavy';
    [mBsc, mExt] = MolModel(height, laserWL, 'meteor', 'standard_atmosphere');
    [aBsc, aExt] = AerModel(height, 'scene', 'marine-moderate');
    [fBsc, fExt] = SeaFogModel(distArr, laserWL, 'scene', seaFogType);
    dataSim = LISMO_Model(distArr, 'tBsc', mBsc + aBsc + fBsc, 'tExt', mExt + aExt + fExt, 'eleAngle', eleAngle, 'laserWL', laserWL, 'visible', 'off');
    save(fullfile(savePath, sprintf('%s.mat', seaFogType)), 'distArr', 'eleAngle', 'aBsc', 'aExt', 'fBsc', 'fExt', 'mBsc', 'mExt', 'dataSim');

end

%% Simulate forward/backward retrievals
load(fullfile(savePath, 'sea-fog-none.mat'));
refH_backward = [9000, 9500];
visibility = ext2vis(mExt * (laserWL / 550)^4 + aExt * (laserWL / 550)^0.5);

% backward retrieval (with well-known refBeta)
refIdx = (distArr >= refH_backward(1)) & (distArr <= refH_backward(2));
[aBsc1, ~] = fernald(distArr, dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR, (dataSim.Nb_FR + dataSim.Nd_FR) * ones(size(distArr)), 20, refH_backward, mean(aBsc(refIdx)), mBsc);
aExt1 = aBsc1 * 20;
visibility1 = ext2vis(aExt1 * (laserWL / 550)^0.5 + mExt * (laserWL / 550)^4);

% backward retrieval (without refBeta)
[aBsc2, ~] = fernald(distArr, dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR, (dataSim.Nb_FR + dataSim.Nd_FR) * ones(size(distArr)), 20, refH_backward, 0, mBsc);
aExt2 = aBsc2 * 20;
visibility2 = ext2vis(aExt2 * (laserWL / 550)^0.5 + mExt * (laserWL / 550)^4);

% forward retrieval
refH_forward = [80, 100];
refHIdx = (distArr >= refH_forward(1)) & (distArr <= refH_forward(2));
sigNoBg = dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR;
lidarConst = mean(sigNoBg(refHIdx)) * (mean(distArr(refHIdx))^2) / mean(aBsc(refHIdx) + mBsc(refHIdx) + fBsc(refHIdx));
[~, aExt3] = quasiRetrieval(distArr, sigNoBg .* distArr'.^2 / lidarConst, mExt', mBsc', 20, 'nIters', 50);
aExt3 = transpose(aExt3);
visibility3 = ext2vis(aExt3 * (laserWL / 550)^0.5 + mExt * (laserWL / 550)^4);

% data visualization
% signal + extinction (with relative error)
figPos = subfigPos([0.08, 0.14, 0.88, 0.83], 1, 3, 0.04, 0);

figure('Position', [0, 10, 600, 400], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', figPos(1, :), 'Units', 'normalized');
rcsTmp = sigNoBg' .* distArr.^2;
rcsTmp(rcsTmp <= 0) = NaN;
semilogx(rcsTmp, distArr / 1e3, '-r'); hold on;

xlim([1e9, 1e12]);
ylim(hRange);

xlabel('Range corrected signal');
ylabel('Distance (km)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02], 'FontSize', 10);

subplot('Position', figPos(2, :), 'Units', 'normalized');

p4 = plot(aExt3 * 1e6, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 1, 'DisplayName', 'Forward');
hold on;

p2 = plot(aExt2 * 1e6, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 1, 'DisplayName', 'Backward (ref. 0)');
p3 = plot(aExt1 * 1e6, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 1, 'DisplayName', 'Backward');
p1 = plot((aExt + fExt) * 1e6, distArr / 1e3, 'color', [53, 46, 51] / 255, 'LineWidth', 1, 'LineStyle', '--', 'DisplayName', 'True Val.');

xlim([0, 200]);
ylim(hRange);

xlabel('Extinction Coeff. (Mm^{-1})');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

legend([p1, p2, p3, p4], 'Location', 'NorthEast');

subplot('Position', figPos(3, :), 'Units', 'normalized');

plot((aExt3 - aExt - fExt) ./ (aExt + fExt) * 100, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 1);
hold on;
plot((aExt1 - aExt - fExt) ./ (aExt + fExt) * 100, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 1);
plot((aExt2 - aExt - fExt) ./ (aExt + fExt) * 100, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 1);

plot([15, 15], [0, 100], '-.k');
plot([-15, -15], [0, 100], '-.k');

plot([0, 0], [0, 100], 'k--');

xlim([-50, 50]);
ylim(hRange);

xlabel('Rel. Err. (%)');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'retrieval-cmp-extinction-no-sea-fog.png'), '-r300');

% signal + visibility (with relative error)
figPos = subfigPos([0.08, 0.14, 0.88, 0.83], 1, 3, 0.04, 0);

figure('Position', [0, 10, 600, 400], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', figPos(1, :), 'Units', 'normalized');
rcsTmp = sigNoBg' .* distArr.^2;
rcsTmp(rcsTmp <= 0) = NaN;
semilogx(rcsTmp, distArr / 1e3, '-r'); hold on;

% xlim();
ylim(hRange);

xlabel('Range corrected signal');
ylabel('Distance (km)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02], 'FontSize', 10);

subplot('Position', figPos(2, :), 'Units', 'normalized');

p4 = semilogx(visibility3, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 1, 'DisplayName', 'Forward');
hold on;

p2 = semilogx(visibility2, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 1, 'DisplayName', 'Backward (ref. 0)');
p3 = semilogx(visibility1, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 1, 'DisplayName', 'Backward');
p1 = semilogx(visibility, distArr / 1e3, 'color', [53, 46, 51] / 255, 'LineWidth', 1, 'LineStyle', '--', 'DisplayName', 'True Val.');

xlim([1, 1e7]);
ylim(hRange);

xlabel('Visibility (m)');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

legend([p1, p2, p3, p4], 'Location', 'NorthEast');

subplot('Position', figPos(3, :), 'Units', 'normalized');

plot((visibility3 - visibility) ./ visibility * 100, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 1);
hold on;
plot((visibility1 - visibility) ./ visibility * 100, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 1);
plot((visibility2 - visibility) ./ visibility * 100, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 1);

plot([0, 0], [0, 100], 'k--');

plot([30, 30], [0, 100], '-.k');
plot([-30, -30], [0, 100], '-.k');

xlim([-50, 50]);
ylim(hRange);

xlabel('Rel. Err. (%)');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'retrieval-cmp-visibility-no-sea-fog.png'), '-r300');

%% Simulate improvements after iterations (sea-fog conditions)
load(fullfile(savePath, 'sea-fog-moderate.mat'));
visibility = ext2vis(mExt * (laserWL / 550)^4 + aExt * (laserWL / 550)^0.5);
nLoops = 50;

% forward retrieval
refH_forward = [80, 100];
lr = (aExt + fExt) ./ (aBsc + fBsc);
refHIdx = (distArr >= refH_forward(1)) & (distArr <= refH_forward(2));
sigNoBg = dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR;
lidarConst = mean(transpose(sigNoBg(refHIdx)) .* distArr(refHIdx).^2) / mean(aBsc(refHIdx) + mBsc(refHIdx) + fBsc(refHIdx)) * exp(2 * mean((aBsc(refHIdx) + fBsc(refHIdx)) .* lr(refHIdx)) * mean(refH_forward));

aBsc4 = NaN(length(distArr), nLoops);
aExt4 = NaN(length(distArr), nLoops);
visibility4 = NaN(length(distArr), nLoops);

for iLoop = 1:nLoops
    [aBscTmp, aExtTmp] = quasiRetrieval(distArr, sigNoBg .* distArr'.^2 / lidarConst, mExt', mBsc', transpose(lr), 'nIters', iLoop);
    visibilityTmp = ext2vis(aExtTmp' * (laserWL / 550)^0.5 + mExt * (laserWL / 550)^4);

    aBsc4(:, iLoop) = aBscTmp;
    aExt4(:, iLoop) = aExtTmp;
    visibility4(:, iLoop) = visibilityTmp;
end

% data visualization
% signal + extinction (with relative error)
figure('Position', [0, 10, 500, 400], 'Units', 'Pixels', 'Color', 'w');

p1 = pcolor(1:nLoops, distArr / 1e3, (aExt4 - repmat(transpose(aExt + fExt), 1, nLoops)) ./ repmat(transpose(aExt + fExt), 1, nLoops) * 100); hold on;
p1.EdgeColor = 'None';

xlim([1, nLoops]);
ylim(hRange);
caxis([-20, 20]);

xlabel('Iteration');
ylabel('Distance (km)');
title('Improvement of Forward Retrieval After Iterations');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02], 'FontSize', 10, 'TickDir', 'out');

cb = colorbar();
titleHandle = get(cb, 'Title');
set(titleHandle, 'String', '[%]');
set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'improvement-error-forward-retrieval-after-iteration.png'), '-r300');

% line plot
colors = colormap('parula');
figure('Position', [0, 10, 500, 400], 'Units', 'Pixels', 'Color', 'w');
figPos = subfigPos([0.08, 0.14, 0.88, 0.83], 1, 2, 0.06, 0);

subplot('Position', figPos(1, :), 'units', 'normalized');
lineInstances = cell(0);

loopArr = [1, 2, 3, 4, 5, 8, 10, 15];
for iLoop = 1:length(loopArr)
    p1 = plot(aExt4(:, loopArr(iLoop)) * 1e6, distArr / 1e3, 'color', colors(floor((iLoop / length(loopArr) * 63) + 1), :), 'DisplayName', sprintf('Iters: %d', loopArr(iLoop))); hold on;
    lineInstances = cat(2, lineInstances, p1);
end

p1 = plot((aExt + fExt) * 1e6, distArr / 1e3, 'color', [53, 46, 51] / 255, 'LineWidth', 1, 'LineStyle', '--', 'DisplayName', 'True Val.'); hold on;
lineInstances = cat(2, lineInstances, p1);

xlim([0, 1.5e3]);
ylim(hRange);

xlabel('Extinction Coeff. (m^{-1})');
ylabel('Height (km)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02], 'FontSize', 10);

legend(lineInstances, 'Location', 'NorthEast');

subplot('Position', figPos(2, :), 'units', 'normalized');
lineInstances = cell(0);

for iLoop = 1:length(loopArr)
    p1 = plot((transpose(aExt4(:, loopArr(iLoop))) - aExt - fExt) ./ (aExt + fExt) * 100, distArr / 1e3, 'color', colors(floor((iLoop / length(loopArr) * 63) + 1), :), 'DisplayName', sprintf('Iters: %d', loopArr(iLoop))); hold on;
    lineInstances = cat(2, lineInstances, p1);
end

plot([15, 15], [0, 100], '-.k');
plot([-15, -15], [0, 100], '-.k');

plot([0, 0], [0, 100], 'k--');

xlim([-50, 50]);
ylim(hRange);

xlabel('Rel. Err. (%)');
ylabel('');

set(gca, 'YTickLabel', '', 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02], 'FontSize', 10);

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'improvement-error-forward-retrieval-after-iteration-line.png'), '-r300');
export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'improvement-error-forward-retrieval-after-iteration-line.pdf'));

%% Simulate forward/backward retrievals under sea-fog (weak-moderate-heavy)
refH_backward = [14500, 15000];

% weak
load(fullfile(savePath, 'sea-fog-weak.mat'));
visibility = ext2vis(mExt * (laserWL / 550)^4 + (aExt + fExt) * (laserWL / 550)^0.5);

lr = (aExt + fExt) ./ (aBsc + fBsc);

% backward retrieval (with well-known refBeta)
refIdx = (distArr >= refH_backward(1)) & (distArr <= refH_backward(2));
[aBsc_Bcw1, ~] = fernald(distArr, dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR, (dataSim.Nb_FR + dataSim.Nd_FR) * ones(size(distArr)), lr, refH_backward, mean(aBsc(refIdx)), mBsc, 5);
aExt_Bsw1 = aBsc_Bcw1 .* lr;
visibility1 = ext2vis(aExt_Bsw1 * (laserWL / 550)^0.5 + mExt * (laserWL / 550)^4);

% backward retrieval (without refBeta)
[aBsc_Bsc2, ~] = fernald(distArr, dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR, (dataSim.Nb_FR + dataSim.Nd_FR) * ones(size(distArr)), lr, refH_backward, 0, mBsc, 5);
aExt_Bsw2 = aBsc_Bsc2 .* lr;
visibility2 = ext2vis(aExt_Bsw2 * (laserWL / 550)^0.5 + mExt * (laserWL / 550)^4);

% forward retrieval
refH_forward = [80, 100];
refHIdx = (distArr >= refH_forward(1)) & (distArr <= refH_forward(2));
sigNoBg = dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR;
lidarConst = mean(sigNoBg(refHIdx)) * (mean(distArr(refHIdx))^2) / mean(aBsc(refHIdx) + mBsc(refHIdx) + fBsc(refHIdx));
[~, aExt3] = quasiRetrieval(distArr, sigNoBg .* distArr'.^2 / lidarConst, mExt', mBsc', lr', 'nIters', 15);
aExt3 = transpose(aExt3);
visibility3 = ext2vis(aExt3 * (laserWL / 550)^0.5 + mExt * (laserWL / 550)^4);

% data visualization
% signal + extinction (with relative error)
figPos = subfigPos([0.08, 0.14, 0.88, 0.83], 1, 3, 0.04, 0);

figure('Position', [0, 10, 600, 400], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', figPos(1, :), 'Units', 'normalized');
rcsTmp = sigNoBg' .* distArr.^2;
rcsTmp(rcsTmp <= 0) = NaN;
semilogx(rcsTmp, distArr / 1e3, '-r'); hold on;

% xlim([1e9, 1e12]);
ylim(hRange);

xlabel('Range corrected signal');
ylabel('Distance (km)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02], 'FontSize', 10);

subplot('Position', figPos(2, :), 'Units', 'normalized');

p4 = plot(aExt3 * 1e6, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 2, 'DisplayName', 'Forward');
hold on;

p2 = plot(aExt_Bsw2 * 1e6, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 2, 'DisplayName', 'Backward (ref. 0)');
% p3 = plot(aExt_Bsw1 * 1e6, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 2, 'DisplayName', 'Backward');
p1 = plot((aExt + fExt) * 1e6, distArr / 1e3, 'color', [53, 46, 51] / 255, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'True Val.');

xlim([0, 3000]);
ylim(hRange);

xlabel('Extinction Coeff. (Mm^{-1})');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

legend([p1, p2, p4], 'Location', 'NorthEast');

subplot('Position', figPos(3, :), 'Units', 'normalized');

plot((aExt3 - aExt - fExt) ./ (aExt + fExt) * 100, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 2);
hold on;
plot((aExt_Bsw2 - aExt - fExt) ./ (aExt + fExt) * 100, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 2);
% p1 = plot((aExt_Bsw1 - aExt - fExt) ./ (aExt + fExt) * 100, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 2);

plot([15, 15], [0, 100], '-.k');
plot([-15, -15], [0, 100], '-.k');

plot([0, 0], [0, 100], 'k--');

xlim([-50, 50]);
ylim(hRange);

xlabel('Rel. Err. (%)');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'retrieval-cmp-extinction-sea-fog-weak.png'), '-r300');

% signal + visibility (with relative error)
figPos = subfigPos([0.08, 0.14, 0.88, 0.83], 1, 3, 0.04, 0);

figure('Position', [0, 10, 600, 400], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', figPos(1, :), 'Units', 'normalized');
rcsTmp = sigNoBg' .* distArr.^2;
rcsTmp(rcsTmp <= 0) = NaN;
semilogx(rcsTmp, distArr / 1e3, '-r'); hold on;

% xlim();
ylim(hRange);

xlabel('Range corrected signal');
ylabel('Distance (km)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02], 'FontSize', 10);

subplot('Position', figPos(2, :), 'Units', 'normalized');

p4 = semilogx(visibility3, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 2, 'DisplayName', 'Forward');
hold on;

p2 = semilogx(visibility2, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 2, 'DisplayName', 'Backward (ref. 0)');
% p3 = semilogx(visibility1, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 2, 'DisplayName', 'Backward');
p1 = semilogx(visibility, distArr / 1e3, 'color', [53, 46, 51] / 255, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'True Val.');

xlim([1, 1e7]);
ylim(hRange);

xlabel('Visibility (m)');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

legend([p1, p2, p4], 'Location', 'NorthEast');

subplot('Position', figPos(3, :), 'Units', 'normalized');

plot((visibility3 - visibility) ./ visibility * 100, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 2);
hold on;
% p2 = plot((visibility1 - visibility) ./ visibility * 100, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 2);
plot((visibility2 - visibility) ./ visibility * 100, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 2);

plot([0, 0], [0, 100], 'k--');

plot([30, 30], [0, 100], '-.k');
plot([-30, -30], [0, 100], '-.k');

xlim([-50, 50]);
ylim(hRange);

xlabel('Rel. Err. (%)');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'retrieval-cmp-visibility-sea-fog-weak.png'), '-r300');

% moderate
load(fullfile(savePath, 'sea-fog-moderate.mat'));
visibility = ext2vis(mExt * (laserWL / 550)^4 + (aExt + fExt) * (laserWL / 550)^0.5);

lr = (aExt + fExt) ./ (aBsc + fBsc);

% backward retrieval (with well-known refBeta)
refIdx = (distArr >= refH_backward(1)) & (distArr <= refH_backward(2));
[aBsc_Bcw1, ~] = fernald(distArr, dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR, (dataSim.Nb_FR + dataSim.Nd_FR) * ones(size(distArr)), lr, refH_backward, mean(aBsc(refIdx)), mBsc, 5);
aExt_Bsw1 = aBsc_Bcw1 .* lr;
visibility1 = ext2vis(aExt_Bsw1 * (laserWL / 550)^0.5 + mExt * (laserWL / 550)^4);

% backward retrieval (without refBeta)
[aBsc_Bsc2, ~] = fernald(distArr, dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR, (dataSim.Nb_FR + dataSim.Nd_FR) * ones(size(distArr)), lr, refH_backward, 0, mBsc, 5);
aExt_Bsw2 = aBsc_Bsc2 .* lr;
visibility2 = ext2vis(aExt_Bsw2 * (laserWL / 550)^0.5 + mExt * (laserWL / 550)^4);

% forward retrieval
refH_forward = [80, 100];
refHIdx = (distArr >= refH_forward(1)) & (distArr <= refH_forward(2));
sigNoBg = dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR;
lidarConst = mean(sigNoBg(refHIdx)) * (mean(distArr(refHIdx))^2) / mean(aBsc(refHIdx) + mBsc(refHIdx) + fBsc(refHIdx));
[~, aExt3] = quasiRetrieval(distArr, sigNoBg .* distArr'.^2 / lidarConst, mExt', mBsc', lr', 'nIters', 50);
aExt3 = transpose(aExt3);
visibility3 = ext2vis(aExt3 * (laserWL / 550)^0.5 + mExt * (laserWL / 550)^4);

% data visualization
% signal + extinction (with relative error)
figPos = subfigPos([0.08, 0.14, 0.88, 0.83], 1, 3, 0.04, 0);

figure('Position', [0, 10, 600, 400], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', figPos(1, :), 'Units', 'normalized');
rcsTmp = sigNoBg' .* distArr.^2;
rcsTmp(rcsTmp <= 0) = NaN;
semilogx(rcsTmp, distArr / 1e3, '-r'); hold on;

% xlim([1e9, 1e12]);
ylim(hRange);

xlabel('Range corrected signal');
ylabel('Distance (km)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02], 'FontSize', 10);

subplot('Position', figPos(2, :), 'Units', 'normalized');

p4 = plot(aExt3 * 1e6, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 2, 'DisplayName', 'Forward');
hold on;

p2 = plot(aExt_Bsw2 * 1e6, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 2, 'DisplayName', 'Backward (ref. 0)');
% p3 = plot(aExt_Bsw1 * 1e6, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 2, 'DisplayName', 'Backward');
p1 = plot((aExt + fExt) * 1e6, distArr / 1e3, 'color', [53, 46, 51] / 255, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'True Val.');

xlim([0, 1500]);
ylim(hRange);

xlabel('Extinction Coeff. (Mm^{-1})');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

legend([p1, p2, p4], 'Location', 'NorthEast');

subplot('Position', figPos(3, :), 'Units', 'normalized');

plot((aExt3 - aExt - fExt) ./ (aExt + fExt) * 100, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 2);
hold on;
plot((aExt_Bsw2 - aExt - fExt) ./ (aExt + fExt) * 100, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 2);
% p1 = plot((aExt_Bsw1 - aExt - fExt) ./ (aExt + fExt) * 100, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 2);

plot([15, 15], [0, 100], '-.k');
plot([-15, -15], [0, 100], '-.k');

plot([0, 0], [0, 100], 'k--');

xlim([-50, 50]);
ylim(hRange);

xlabel('Rel. Err. (%)');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'retrieval-cmp-extinction-sea-fog-moderate.png'), '-r300');

% signal + visibility (with relative error)
figPos = subfigPos([0.08, 0.14, 0.88, 0.83], 1, 3, 0.04, 0);

figure('Position', [0, 10, 600, 400], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', figPos(1, :), 'Units', 'normalized');
rcsTmp = sigNoBg' .* distArr.^2;
rcsTmp(rcsTmp <= 0) = NaN;
semilogx(rcsTmp, distArr / 1e3, '-r'); hold on;

% xlim();
ylim(hRange);

xlabel('Range corrected signal');
ylabel('Distance (km)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02], 'FontSize', 10);

subplot('Position', figPos(2, :), 'Units', 'normalized');

p4 = semilogx(visibility3, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 2, 'DisplayName', 'Forward');
hold on;

p2 = semilogx(visibility2, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 2, 'DisplayName', 'Backward (ref. 0)');
% p3 = semilogx(visibility1, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 2, 'DisplayName', 'Backward');
p1 = semilogx(visibility, distArr / 1e3, 'color', [53, 46, 51] / 255, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'True Val.');

xlim([1, 1e7]);
ylim(hRange);

xlabel('Visibility (m)');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

legend([p1, p2, p4], 'Location', 'NorthEast');

subplot('Position', figPos(3, :), 'Units', 'normalized');

hold on;
plot((visibility3 - visibility) ./ visibility * 100, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 2);
% p2 = plot((visibility1 - visibility) ./ visibility * 100, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 2);
plot((visibility2 - visibility) ./ visibility * 100, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 2);

plot([0, 0], [0, 100], 'k--');

plot([30, 30], [0, 100], '-.k');
plot([-30, -30], [0, 100], '-.k');
hold off;

xlim([-50, 50]);
ylim(hRange);

xlabel('Rel. Err. (%)');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'retrieval-cmp-visibility-sea-fog-moderate.png'), '-r300');

% heavy
load(fullfile(savePath, 'sea-fog-heavy.mat'));
visibility = ext2vis(mExt * (laserWL / 550)^4 + (aExt + fExt) * (laserWL / 550)^0.5);

lr = (aExt + fExt) ./ (aBsc + fBsc);

% backward retrieval (with well-known refBeta)
refIdx = (distArr >= refH_backward(1)) & (distArr <= refH_backward(2));
[aBsc_Bcw1, ~] = fernald(distArr, dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR, (dataSim.Nb_FR + dataSim.Nd_FR) * ones(size(distArr)), lr, refH_backward, mean(aBsc(refIdx)), mBsc, 5);
aExt_Bsw1 = aBsc_Bcw1 .* lr;
visibility1 = ext2vis(aExt_Bsw1 * (laserWL / 550)^0.5 + mExt * (laserWL / 550)^4);

% backward retrieval (without refBeta)
[aBsc_Bsc2, ~] = fernald(distArr, dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR, (dataSim.Nb_FR + dataSim.Nd_FR) * ones(size(distArr)), lr, refH_backward, 0, mBsc, 5);
aExt_Bsw2 = aBsc_Bsc2 .* lr;
visibility2 = ext2vis(aExt_Bsw2 * (laserWL / 550)^0.5 + mExt * (laserWL / 550)^4);

% forward retrieval
refH_forward = [80, 100];
refHIdx = (distArr >= refH_forward(1)) & (distArr <= refH_forward(2));
sigNoBg = dataSim.N_FR_Poiss - dataSim.Nb_FR - dataSim.Nd_FR;
lidarConst = mean(sigNoBg(refHIdx)) * (mean(distArr(refHIdx))^2) / mean(aBsc(refHIdx) + mBsc(refHIdx) + fBsc(refHIdx));
[~, aExt3] = quasiRetrieval(distArr, sigNoBg .* distArr'.^2 / lidarConst, mExt', mBsc', lr', 'nIters', 100);
aExt3 = transpose(aExt3);
visibility3 = ext2vis(aExt3 * (laserWL / 550)^0.5 + mExt * (laserWL / 550)^4);

% data visualization
% signal + extinction (with relative error)
figPos = subfigPos([0.08, 0.14, 0.88, 0.83], 1, 3, 0.04, 0);

figure('Position', [0, 10, 600, 400], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', figPos(1, :), 'Units', 'normalized');
rcsTmp = sigNoBg' .* distArr.^2;
rcsTmp(rcsTmp <= 0) = NaN;
semilogx(rcsTmp, distArr / 1e3, '-r'); hold on;

% xlim([]);
ylim(hRange);

xlabel('Range corrected signal');
ylabel('Distance (km)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02], 'FontSize', 10);

subplot('Position', figPos(2, :), 'Units', 'normalized');

p4 = plot(aExt3 * 1e6, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 2, 'DisplayName', 'Forward');
hold on;

p2 = plot(aExt_Bsw2 * 1e6, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 2, 'DisplayName', 'Backward (ref. 0)');
% p3 = plot(aExt_Bsw1 * 1e6, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 2, 'DisplayName', 'Backward');
p1 = plot((aExt + fExt) * 1e6, distArr / 1e3, 'color', [53, 46, 51] / 255, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'True Val.');

xlim([0, 3000]);
ylim(hRange);

xlabel('Extinction Coeff. (Mm^{-1})');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

legend([p1, p2, p4], 'Location', 'NorthEast');

subplot('Position', figPos(3, :), 'Units', 'normalized');

plot((aExt3 - aExt - fExt) ./ (aExt + fExt) * 100, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 2);
hold on;
plot((aExt_Bsw2 - aExt - fExt) ./ (aExt + fExt) * 100, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 2);
% p1 = plot((aExt_Bsw1 - aExt - fExt) ./ (aExt + fExt) * 100, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 2);

plot([15, 15], [0, 100], '-.k');
plot([-15, -15], [0, 100], '-.k');

plot([0, 0], [0, 100], 'k--');

xlim([-50, 50]);
ylim(hRange);

xlabel('Rel. Err. (%)');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'retrieval-cmp-extinction-sea-fog-heavy.png'), '-r300');

% signal + visibility (with relative error)
figPos = subfigPos([0.08, 0.14, 0.88, 0.83], 1, 3, 0.04, 0);

figure('Position', [0, 10, 600, 400], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', figPos(1, :), 'Units', 'normalized');
rcsTmp = sigNoBg' .* distArr.^2;
rcsTmp(rcsTmp <= 0) = NaN;
semilogx(rcsTmp, distArr / 1e3, '-r'); hold on;

% xlim();
ylim(hRange);

xlabel('Range corrected signal');
ylabel('Distance (km)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02], 'FontSize', 10);

subplot('Position', figPos(2, :), 'Units', 'normalized');

p4 = semilogx(visibility3, distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 2, 'DisplayName', 'Forward');
hold on;

p2 = semilogx(visibility2, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 2, 'DisplayName', 'Backward (ref. 0)');
% p3 = semilogx(visibility1, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 2, 'DisplayName', 'Backward');
p1 = semilogx(visibility, distArr / 1e3, 'color', [53, 46, 51] / 255, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'True Val.');

xlim([1, 1e7]);
ylim(hRange);

xlabel('Visibility (m)');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

legend([p1, p2, p4], 'Location', 'NorthEast');

subplot('Position', figPos(3, :), 'Units', 'normalized');

p3 = plot(smooth((visibility3 - visibility) ./ visibility * 100, 8), distArr / 1e3, 'color', [238, 107, 171] / 255, 'LineWidth', 2);
hold on;
% p2 = plot((visibility1 - visibility) ./ visibility * 100, distArr / 1e3, 'color', [20, 110, 178] / 255, 'LineWidth', 2);
p1 = plot((visibility2 - visibility) ./ visibility * 100, distArr / 1e3, 'color', [201, 229, 240] / 255, 'LineWidth', 2);

plot([0, 0], [0, 100], 'k--');

plot([30, 30], [0, 100], '-.k');
plot([-30, -30], [0, 100], '-.k');

xlim([-50, 50]);
ylim(hRange);

xlabel('Rel. Err. (%)');
ylabel('');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'YTickLabel', '', 'TickLen', [0.02, 0.02], 'FontSize', 10);

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', 'retrieval-cmp-visibility-sea-fog-heavy.png'), '-r300');