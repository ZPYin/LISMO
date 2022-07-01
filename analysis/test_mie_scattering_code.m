clc;
close all;

%% Parameter definition
wl = 200:100:12000;   % wavelenght of incidient light. (nm)
radius = 8;   % radius of water droplet. (micron)

%% Refrective index of water
refIdxReal = zeros(size(wl));
refIdxImg = zeros(size(wl));
for iWL = 1:length(wl)
    [thisRefIdxReal, thisRefIdxImg] = refIdxWater(wl(iWL));

    refIdxReal(iWL) = thisRefIdxReal;
    refIdxImg(iWL) = thisRefIdxImg;
end

%% Mie scattering
qe = zeros(size(wl));   % efficiencies for extinction
qs = zeros(size(wl));   % efficiencies for scattering
qa = zeros(size(wl));   % efficiencies for absorption

for iWL = 1:length(wl)
    m = refIdxReal(iWL) + 1i * refIdxImg(iWL);
    k0 = 2 * pi / wl(iWL);
    a = radius * 1000;
    res = Mie(m, k0 * a);

    qe(iWL) = res(1);
    qs(iWL) = res(2);
    qa(iWL) = res(3);
end

%% Data visualization
figure('Position', [0, 10, 550, 300], 'Units', 'Pixels', 'Color', 'w');

p1 = plot(wl/1e3, qe, '-k', 'DisplayName', 'Q_e'); hold on;
p2 = plot(wl/1e3, qs, '-r', 'DisplayName', 'Q_s');
p3 = plot(wl/1e3, qa, '-g', 'DisplayName', 'Q_a');

xlabel('\lambda (\mum)');
ylabel('Efficiency');

xlim([0, 12]);
ylim([0, 3.5]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.01]);

legend([p1, p2, p3], 'Location', 'NorthEast');

export_fig(gcf, 'D:\Coding\Matlab\LISMO\image\efficiency_water.png', '-r300');