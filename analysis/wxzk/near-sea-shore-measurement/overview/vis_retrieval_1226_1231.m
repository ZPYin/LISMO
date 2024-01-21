% ��������ʱ�䷶Χ�ڵİ���ɨ���ܼ����״����ݣ����������ǰ��ɢ���ܼ����ǽ��жԱ�
%
% ����ƽ
% 2024-01-21

clc;
close all;

%% Parameter Definition
dataFolder = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\RawData\QingDao1';   % ������Ŀ¼
tRange = [datenum(2023, 12, 25, 0, 0, 0), datenum(2023, 12, 31, 23, 59, 59)];   % ��������ʱ�䷶Χ
distOffset = -48.75;   % ���У��
iCh = 1;   % ��ȡͨ���±꣨һ���ĸ�ͨ��������������Ϊ������Զ��ͨ����
flagAve4Vis = false;   % �Ƿ��ÿ��ɨ�����ڽ��������ۼ�
iPrfInScan4Vis = 1;   % ��������������ۼӣ���ȡÿ��ɨ�������еĵ�iPrfInScan4Vis������
flagOverlapCor = true;   % �Ƿ�����ص�����У��
olHeight = 1000;   % ä���߶�(��)����Ҫ�����������㷨
visRetMethod = 'xian';   % �ܼ��ȷ��ݷ���������ʹ�������鷽��
overlapFile = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\wxzk\near-sea-shore-measurement\2023-12-29\overlap_20231229.mat';   % �ص����������ļ����ص����ӿ���ͨ��overlap_estimation_xxxx.m������й���
visFile = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\analysis\wxzk\near-sea-shore-measurement\vis-fixed-platform.mat';   % �ܼ��������ļ�

%% Read Data
fullData = struct();
fullData.mTime = [];
fullData.rawSignal = [];
fullData.nShots = [];
for thisDay = floor(tRange(1)):floor(tRange(2))
    dayPath = fullfile(dataFolder, datestr(thisDay, 'yyyy'), datestr(thisDay, 'mm'), datestr(thisDay, 'dd'));
    scanPaths = listdir(dayPath, '.*', 1);

    for iScan = 1:length(scanPaths)
        fprintf('Finished %6.2f%%: reading %s\n', (iScan - 1) / length(scanPaths) * 100, scanPaths{iScan});

        scanFiles = listfile(scanPaths{iScan}, '\w*', 1);
        if ~ flagAve4Vis
            scanFiles = scanFiles{iPrfInScan4Vis};
        end

        lData = readVIS(scanFiles, 'isDir', false);

        fullData.mTime = cat(2, fullData.mTime, mean(lData.startTime));
        fullData.nShots = cat(2, fullData.nShots, sum(lData.nShots));
        fullData.rawSignal = cat(1, fullData.rawSignal, sum(lData.rawSignal, 1));
    end
end

% Read Vis Data
visData = load(visFile);

%% Preprocess
range = ((1:lData.nBins(1)) - 0.5) * lData.hRes(1) + distOffset;

fullData.bg = nanmean(squeeze(fullData.rawSignal(:, iCh, (end - 20):end)), 2);
fullData.signal = squeeze(fullData.rawSignal(:, iCh, :)) - repmat(fullData.bg, 1, lData.nBins(1));
fullData.rcs = fullData.signal .* repmat(range, size(fullData.rawSignal, 1), 1).^2;
fullData.snr = fullData.signal ./ sqrt(squeeze(fullData.rawSignal(:, iCh, :)));
fullData.lowSNRMask = false(size(fullData.signal));
for iPrf = 1:size(fullData.signal, 1)
    snrPrf = fullData.snr(iPrf, :);
    snrPrf(range < olHeight) = NaN;
    snrLow = find(snrPrf < 1, 1);
    fullData.lowSNRMask(iPrf, snrLow:end) = true;
end

%% Overlap Correction
if flagOverlapCor
    ol = load(overlapFile);

    if exist('fullData', 'var')
        fullData.signal = fullData.signal ./ repmat(transpose(ol.ov), size(fullData.signal, 1), 1);
        fullData.rcs = fullData.rcs ./ repmat(transpose(ol.ov), size(fullData.rcs, 1), 1);
    end

    if exist('scanData', 'var')
        scanData.signal = scanData.signal ./ repmat(transpose(ol.ov), size(scanData.signal, 1), 1);
        scanData.rcs = scanData.rcs ./ repmat(transpose(ol.ov), size(scanData.rcs, 1), 1);
    end
end

%% extinction&visibility retrieval
fullData.ext = NaN(size(fullData.signal));
fullData.vis = NaN(size(fullData.signal));

for iPrf = 1:size(fullData.signal, 1)

    fullData.ext(iPrf, :) = extRet_Xian(range, fullData.signal(iPrf, :), fullData.bg(iPrf), 'minSNR', 0.5, 'rangeFullOverlap', 200);
end

fullData.ext(fullData.lowSNRMask) = NaN;
fullData.vis = ext2vis(fullData.ext);
fullData.vis(isnan(fullData.vis)) = 1e5;

%% Display

% full data rcs
figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
hold on;
rcsTmp = fullData.rcs;
rcsTmp(fullData.lowSNRMask) = NaN;
p1 = pcolor(fullData.mTime, range * 1e-3, transpose(rcsTmp));
p1.EdgeColor = 'none';
caxis([0, 0.1e11]);
colormap('jet');

xlabel('ʱ��');
ylabel('���� (ǧ��)');
title(sprintf('��������״�ʱ���ź�ͼ (%s��%s)', datestr(fullData.mTime(1), 'yyyy-mm-dd'), datestr(fullData.mTime(end), 'mm-dd')));

xlim([fullData.mTime(1), fullData.mTime(end)]);
ylim([0, 4]);

set(gca, 'XMinorTick', 'on', 'xtick', floor(tRange(1)):ceil(tRange(2)), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 11);
ax = gca;
ax.XAxis.MinorTickValues = floor(tRange(1)):datenum(0, 1, 0, 6, 0, 0):ceil(tRange(2));
datetick(gca, 'x', 'mm-dd', 'keepticks', 'keeplimits');

cb = colorbar('Position', [0.89, 0.2, 0.03, 0.6], 'Units', 'Normalized');
titleHandle = get(cb, 'Title');
set(titleHandle, 'String', '[Mm-1sr-1]');
set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

% full data ext
figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
hold on;
extTmp = fullData.ext;
extTmp(fullData.lowSNRMask) = NaN;
p1 = pcolor(fullData.mTime, range * 1e-3, transpose(extTmp) * 1e3);
p1.EdgeColor = 'none';
caxis([0, 3]);
colormap('jet');

xlabel('ʱ��');
ylabel('���� (ǧ��)');
title(sprintf('����ϵ�� (%s��%s)', datestr(fullData.mTime(1), 'yyyy-mm-dd'), datestr(fullData.mTime(end), 'mm-dd')));

xlim([fullData.mTime(1), fullData.mTime(end)]);
ylim([0, 4]);

set(gca, 'XMinorTick', 'on', 'xtick', floor(tRange(1)):ceil(tRange(2)), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 11);
ax = gca;
ax.XAxis.MinorTickValues = floor(tRange(1)):datenum(0, 1, 0, 6, 0, 0):ceil(tRange(2));
datetick(gca, 'x', 'mm-dd', 'keepticks', 'keeplimits');

cb = colorbar('Position', [0.89, 0.2, 0.03, 0.6], 'Units', 'Normalized');
titleHandle = get(cb, 'Title');
set(titleHandle, 'String', '[km-1]');
set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

% full data vis
figure('Position', [0, 30, 570, 300], 'Units', 'Pixels', 'Color', 'w');

subplot('Position', [0.1, 0.2, 0.75, 0.7], 'Units', 'Normalized');
hold on;
visTmp = fullData.vis;
visTmp(fullData.lowSNRMask) = NaN;
p1 = pcolor(fullData.mTime, range * 1e-3, transpose(visTmp) * 1e-3);
p1.EdgeColor = 'none';
caxis([0, 30]);
load('vis_colormap.mat');
colormap(double(visColorbar) / 255);

xlabel('ʱ��');
ylabel('���� (ǧ��)');
title(sprintf('�ܼ��� (%s��%s)', datestr(fullData.mTime(1), 'yyyy-mm-dd'), datestr(fullData.mTime(end), 'mm-dd')));

xlim([fullData.mTime(1), fullData.mTime(end)]);
ylim([0, 4]);

set(gca, 'XMinorTick', 'on', 'xtick', floor(tRange(1)):ceil(tRange(2)), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'tickdir', 'out', 'fontsize', 11);
ax = gca;
ax.XAxis.MinorTickValues = floor(tRange(1)):datenum(0, 1, 0, 6, 0, 0):ceil(tRange(2));
datetick(gca, 'x', 'mm-dd', 'keepticks', 'keeplimits');

cb = colorbar('Position', [0.89, 0.2, 0.03, 0.6], 'Units', 'Normalized');
titleHandle = get(cb, 'Title');
set(titleHandle, 'String', '[km]');
set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

% full data comparison
figure('Position', [0, 10, 550, 300], 'Units', 'Pixels', 'Color', 'w');

hold on;
s1 = plot(visData.mTime, visData.vis * 1e-3, 'color', 'k', 'Marker', '.', 'MarkerFaceColor', 'k', 'markeredgecolor', 'k', 'DisplayName', '�ܼ�����');
s2 = plot(fullData.mTime, fullData.vis(:, 80) * 1e-3, 'color', 'b', 'Marker', '.', 'MarkerFaceColor', 'b', 'markeredgecolor', 'b', 'DisplayName', sprintf('��������״�(%3d��)', floor(range(80))));
hold off;

xlim([fullData.mTime(1), fullData.mTime(end)]);
ylim([0, 50]);

xlabel('ʱ��');
ylabel('�ܼ��� (ǧ��)');
title(sprintf('�ܼ������жԱ�(%s��%s)', datestr(fullData.mTime(1), 'yyyy-mm-dd'), datestr(fullData.mTime(end), 'mm-dd')));

set(gca, 'XMinorTick', 'on', 'xtick', floor(tRange(1)):ceil(tRange(2)), 'YMinorTick', 'on', 'layer', 'top', 'box', 'on', 'fontsize', 11);
ax = gca;
ax.XAxis.MinorTickValues = floor(tRange(1)):datenum(0, 1, 0, 6, 0, 0):ceil(tRange(2));
datetick(gca, 'x', 'mm-dd', 'keepticks', 'keeplimits');
legend([s1, s2], 'location', 'NorthEast');