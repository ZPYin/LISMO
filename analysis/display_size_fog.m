clc;
close all;

%% Parameter definition
wc = 0.01:0.01:1;
r = 0.001:0.01:45;   % um

%% Size dist function
funAdv = cell(1, length(wc));
funRad = cell(1, length(wc));
for iWC = 1:length(wc)
    funAdv{iWC} = fogSD_KM(wc(iWC), 'tsype', 'advection');
    funRad{iWC} = fogSD_KM(wc(iWC), 'tsype', 'radiation');
end

%% Data visualization
% advection
sd = zeros(length(wc), length(r));
for iWC = 1:length(wc)
    sd(iWC, :) = funAdv{iWC}(r);
end
figure('Position', [0, 10, 400, 300], 'Units', 'Pixels', 'Color', 'w');

p1 = pcolor(r, wc, sd);
p1.EdgeColor = 'none';
colormap('jet');

xlabel('radius (\mum)');
ylabel('Water content (g*m^{-3})');
title('Advection fog');

xlim([0, 45]);
ylim([0, 1]);
caxis([0, 5e6]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'Layer', 'top', 'tickdir', 'out');

colorbar;

export_fig(gcf, 'D:\Coding\Matlab\LISMO\image\advection_fog_size_distribution.png', '-r300');

% radiation
sd = zeros(length(wc), length(r));
for iWC = 1:length(wc)
    sd(iWC, :) = funRad{iWC}(r);
end
figure('Position', [0, 10, 400, 300], 'Units', 'Pixels', 'Color', 'w');

p1 = pcolor(r, wc, sd);
p1.EdgeColor = 'none';
colormap('jet');

xlabel('radius (\mum)');
ylabel('Water content (g*m^{-3})');
title('Radiation fog');

xlim([0, 45]);
ylim([0, 1]);
caxis([0, 1e8]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'Layer', 'top', 'tickdir', 'out');

colorbar;

export_fig(gcf, 'D:\Coding\Matlab\LISMO\image\radiation_fog_size_distribution.png', '-r300');