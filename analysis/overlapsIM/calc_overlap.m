clc;

%% Parameter Definition
d = 0.2;   % distance between laser beam and telescope
r0 = 0.01;   % initial laser width.
theta1 = 0.15e-3;   % divergence of laser beam.
rOF = 2e-4;   % core width of optical fiber
fImg = [0.500, 0.50245, 0.50498];
ox = [-0.00008, -0.001, -0.002];
% fImg = [0.500, 0.500 + (-0.02117 + 0.022), 0.500 + (0.022 - 0.02034), 0.5 + (0.022 - 19.51), 0.5 + (0.022 - 0.01868), 0.5 + (0.022 - 0.01785), 0.5 + (0.022 - 0.01702)];
% ox = [-0.00008, -0.00008 - (0.00071 - 0.00037), -0.00008 - (0.00105 - 0.00037), -0.00008 - (0.00139 - 0.00037), -0.00008 - (0.00173 - 0.00037), -0.00008 - (0.00207 - 0.00037), -0.00008 - (0.00241 - 0.00037)];
% fImg = [0.500, 0.500 + (0.002 - 0.00042), 0.500 + (0.002 + 0.00038), 0.500 + (0.002 + 0.00272)];
% ox = [-0.0002, -0.00160, -0.00253, -0.00513];


%% Overlap calculation
h = logspace(0, log10(1500), 200);
Ov = NaN(length(ox), length(h));

for iOx = 1:length(ox)
    
    for iH = 1:length(h)
        r1 = r0 + h(iH) * theta1 / 2;
        x0 = -d * fImg(iOx) / h(iH) - ox(iOx);
        r2 = r1 * (fImg(iOx) / h(iH));

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

p = plot(h, sum(Ov, 1), '--k', 'DisplayName', 'Overall');
lineInstance = cat(2, lineInstance, p);
plot([0, 1e6], [1, 1], '-.k');
hold off;

xlabel('Distance (m)');
ylabel('Geometrical Overlap Factor [a.u.]');

grid off;

xlim([0, 1500]);
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


export_fig(gcf, sprintf('overlap_function_sim_%dstages.png', length(ox)), '-r300');