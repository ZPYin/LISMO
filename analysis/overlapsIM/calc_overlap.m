clc;

%% Parameter Definition
d = 0.2;   % distance between laser beam and telescope. (m)
r0 = 0.005;   % initial laser radius. (m)
theta1 = 0.05e-3;   % divergence of laser beam. (rad)
rOF = 75e-6;   % core width of optical fiber. (200 um for 3-stage; 150 um for 7-stage)

% 3-stage
% fImg = [0.500, 0.50245, 0.50498];
% ox = [-0.00008, -0.001, -0.002];

% 7-stage
% fImg = [0.500, 0.500 + (-0.02117 + 0.022), 0.500 + (0.022 - 0.02034), 0.5 + (0.022 - 0.01951), 0.5 + (0.022 - 0.01868), 0.5 + (0.022 - 0.01785), 0.5 + (0.022 - 0.01702)];
% ox = [-0.00008, -0.00008 - (0.00071 - 0.00037), -0.00008 - (0.00105 - 0.00037), -0.00008 - (0.00139 - 0.00037), -0.00008 - (0.00173 - 0.00037), -0.00008 - (0.00207 - 0.00037), -0.00008 - (0.00241 - 0.00037)];

% 4-stage
fImg = [0.500, 0.500 + (0.0018 - 0.001634), 0.500 + (0.0018 - 0.001468), 0.500 + (0.0018 - 0.001302)];
ox = [-0.00008, -0.00008 - 0.33e-3, -0.00008 - 0.33e-3 * 2, -0.00008 - 0.33e-3 * 3];

% 13-stage
% fImg = ((1:35) - 1) * 0.366e-3 + 0.5;
% ox = -((1:35) - 1) * 0.15e-3 - 0.00008;

% Display
xRange = [0, 3000];

%% Overlap calculation
h = logspace(0, log10(xRange(2)), 700);
Ov = NaN(length(ox), length(h));

for iOx = 1:length(ox)
    
    for iH = 1:length(h)
        r1 = r0 + h(iH) * theta1 / 2;    % radius of laser beam at height h
        x0 = -d * fImg(iOx) / h(iH) - ox(iOx);   % position of laser beam after projection by the telescope
        r2 = r1 * (fImg(iOx) / h(iH));   % radius of imaged circle after projection

        func = @(r, phi) (1 / (2 * pi * r2^2)) * exp(-0.5 * ((r.*cos(phi) - x0).^2 ./ r2^2 + (r .* sin(phi) / r2).^2)) .* r;
        Ov(iOx, iH) = integral2(func, 0, rOF, 0, 2*pi);
    end

end

%% Display
figure('Position', [0, 10, 550, 350], 'Units', 'Pixels', 'Color', 'w');
hold on;
lineInstance = [];
cmap = colormap('parula');
p = plot(h, Ov(1, :), 'color', 'm', 'LineWidth', 2, 'DisplayName', sprintf('Stage %d', 1));
lineInstance = cat(2, lineInstance, p);
for iStage = 2:length(ox)
    p = plot(h, Ov(iStage, :), 'color', cmap(round(iStage / (length(ox) + 5) * 64), :), 'LineWidth', 2, 'DisplayName', sprintf('Stage %d', iStage));
    lineInstance = cat(2, lineInstance, p);
end

p = plot(h, sum(Ov, 1), '--k', 'LineWidth', 2, 'DisplayName', 'Overall');
lineInstance = cat(2, lineInstance, p);
plot([0, 1e6], [1, 1], '-.k');
hold off;

xlabel('Distance (m)');
ylabel('Geometrical Overlap Factor [a.u.]');

grid off;

xlim(xRange);
ylim([-0.1, 1.2]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 12, 'LineWidth', 2);

legend(lineInstance, 'Location', 'best');

pos = get(gca, 'Position');
axes('Parent', gcf, 'Position', [.4 .3 .25 .25]);
hold on;
plot(h, Ov(1, :), 'color', 'm', 'LineWidth', 2, 'DisplayName', sprintf('Stage %d', 1));
for iStage = 2:length(ox)
    plot(h, Ov(iStage, :), 'color', cmap(round(iStage / (length(ox) + 5) * 64), :), 'LineWidth', 2, 'DisplayName', sprintf('Stage %d', iStage));
end

plot(h, sum(Ov, 1), '--k', 'DisplayName', 'Overall');
plot([0, 1e6], [1, 1], '-.k');
hold off;

xlabel('');
ylabel('');

grid off;

xlim([0, 200]);
ylim([-0.1, 1.2]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'YTickLabel', '', 'LineWidth', 2);

export_fig(gcf, sprintf('overlap_function_sim_%dstages-stage-far-range.pdf', length(ox)));