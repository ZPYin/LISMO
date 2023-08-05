clc;
global LISMO_VARS;

%% ��������
dataFolder = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\data\1030\2023-07-16';
deadtime = 22;
firstRangeBin = 15;
tRange = [datenum(2023, 7, 17, 0, 0, 0), datenum(2023, 7, 17, 0, 30, 0)];

%% ��ȡ����
data = readALADats(dataFolder, 'tRange', tRange);

%% ����Ԥ����
rawSigPCR = data.rawSignal / (50 * data.nShots(1) * 1e-3);
corSigPCR = rawSigPCR ./ (1 - deadtime * rawSigPCR * 1e-3);
corSigPC = corSigPCR * 50 * data.nShots(1) * 1e-3;
bg = nanmean(corSigPC(:, (end - 80):(end - 10), :), 2);
sigPC = corSigPC - repmat(bg, 1, size(corSigPC, 2), 1);
height = ((1:size(sigPC, 2)) - firstRangeBin + 0.5) * data.hRes(1);
rcs = sigPC .* repmat(reshape(height, 1, length(height), 1), size(sigPC, 1), 1, size(sigPC, 3)).^2;

%% ���ݿ��ӻ�

% �����۲�ʱ�ո߶�ͼ
figure('Position', [0, 0, 500, 250], 'Units', 'Pixels', 'Color', 'w');

p1 = pcolor(data.mTime, height / 1e3, squeeze(rcs(2, :, :) / 1e9));
p1.EdgeColor = 'None';

xlabel('ʱ�� (Сʱ:����)');
ylabel('���� (ǧ��)');
title(sprintf('��������״���������ź� %s', datestr(data.mTime(1), 'yyyy-mm-dd')));

caxis([0, 10]);
colormap('jet');

set(gca, 'TickDir', 'out', 'YMinorTick', 'on', 'Box', 'on');
xlim([min(data.mTime), max(data.mTime)]);
ylim([0, 20]);
datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

colorbar();

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', '��Сä������ʹ���ź�����չʾ.png'), '-r300');

%% ��Զ̽�����չʾͼ
figure('Position', [0, 10, 500, 270], 'Units', 'Pixels', 'Color', 'w');

rcsFR = squeeze(nansum(sigPC(2, :, :), 3));
rcsNR = squeeze(nansum(sigPC(1, :, :), 3));

p1 = semilogy(height, rcsFR, '-b', 'LineWidth', 2, 'DisplayName', 'Զ��ͨ��');  hold on;
p2 = semilogy(height, rcsNR, '-g', 'LineWidth', 2, 'DisplayName', '����ͨ��');  hold on;

xlabel('���� (��)');
ylabel('̽���ź� (������)');

xlim([0, 1500]);
ylim([1, 1e7]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.03, 0.02]);

figure('Position', [0, 10, 500, 270], 'Units', 'Pixels', 'Color', 'w');

rcsFR = squeeze(nansum(rcs(2, :, :), 3));
rcsNR = squeeze(nansum(rcs(1, :, :), 3));

p1 = semilogy(height, rcsFR, '-b', 'LineWidth', 2, 'DisplayName', 'Զ��ͨ��');  hold on;
p2 = semilogy(height, rcsNR, '-g', 'LineWidth', 2, 'DisplayName', '����ͨ��');  hold on;

xlabel('���� (��)');
ylabel('���������ź�');

xlim([0, 1500]);
ylim([1e5, 1e10]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.03, 0.02]);
