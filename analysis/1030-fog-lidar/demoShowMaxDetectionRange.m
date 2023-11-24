clc;
global LISMO_VARS;

%% ��������
% dataFolder = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\data\1030\2023-08-02';
dataFolder = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\data\1030\2023-11-07';
deadtime = 22;
firstRangeBin = 16;
% tRange = [datenum(2023, 2, 2, 20, 0, 0), datenum(2023, 8, 3, 5, 0, 0)];
tRange = [datenum(2023, 11, 7, 11, 18, 0), datenum(2023, 11, 7, 11, 20, 0)];

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

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', '�����۲��ź�ʱ�ո߶�ͼ.png'), '-r300');

% ��Զ̽�����չʾͼ
figure('Position', [0, 10, 500, 270], 'Units', 'Pixels', 'Color', 'w');

% tForDetectionRange = [datenum(2023, 7, 18, 20, 10, 0), datenum(2023, 7, 18, 21, 10, 0)];
tForDetectionRange = [datenum(2023, 11, 7, 10, 0, 0), datenum(2023, 11, 7, 13, 30, 0)];
isInAVGRange = (data.mTime >= tForDetectionRange(1)) & (data.mTime <= tForDetectionRange(2));
snr = smooth(sum(sigPC(2, :, isInAVGRange), 3), 4) ./ sqrt(smooth(sum(corSigPC(2, :, isInAVGRange), 3), 4));
snr_nr = smooth(sum(sigPC(1, :, isInAVGRange), 3), 4) ./ sqrt(smooth(sum(corSigPC(1, :, isInAVGRange), 3), 4));
startSearchIdx = 200;
idxSNRlt3 = find(snr(startSearchIdx:end) < 3, 1, 'first') + startSearchIdx;

snr(snr <= 0) = NaN;
p1 = semilogy(height / 1e3, snr, '-b', 'LineWidth', 2);  hold on;
p4 = semilogy(height / 1e3, snr_nr, '-g', 'LineWidth', 2); hold on;
p2 = semilogy([0, 100], [3, 3], '--r');
p3 = scatter(height(idxSNRlt3) / 1e3, snr(idxSNRlt3), 5, 'Marker', 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');


text(height(idxSNRlt3) / 1e3 + 0.5, snr(idxSNRlt3) * 3, sprintf('����: %5.2fǧ��', height(idxSNRlt3) / 1e3), 'Units', 'Data', 'color', 'r', 'fontsize', 15);

xlabel('���� (ǧ��)');
ylabel('�����');
title(sprintf('��������״���Զ̽�����չʾ (%s)', datestr(mean(data.mTime(isInAVGRange)), 'yyyy��mm��dd�� HH:MM')));

xlim([0, 30]);
ylim([1e-1, 1e4]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.03, 0.02]);

export_fig(gcf, fullfile(LISMO_VARS.projectDir, 'image', '��Զ̽�����.png'), '-r300');