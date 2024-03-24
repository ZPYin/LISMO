clc;
close all;
global LISMO_VARS;

%% ��������
lidarPath = 'C:\Users\ZPYin\Documents\Data\1030-fog-lidar\data\1030\2023-11-07';
tRange = [datenum(2023, 11, 7, 11, 40, 0), datenum(2023, 11, 7, 12, 50, 0)];
deadTime = [20, 20];   % ��ʱ�䣬��λns
hOffset = -225;
zenithAngle = 0;
dbFile = '';
cloudMinExt = 200e-6;   % m-1
overlapFile = 'C:\Users\ZPYin\Documents\Coding\Matlab\LISMO\data\1030-fog-lidar-overlap-2023-11-07.mat';
idxVisCmp = 387;   % ���ڶԱȵ��ܼ��Ⱦ������±�
ratioNear = 76;

%% ���ݶ�ȡ
lData = readALADats(lidarPath, 'tRange', tRange, 'nmaxbin', 2000);
ol = load(overlapFile);

%% ����Ԥ����

% ��ʱ��У��
rawSig = lData.rawSignal;
pc2cr = 1 / (1500 / lData.hRes(1) * lData.nShots(1));
rawSigCR = rawSig * pc2cr;
rawSigCRCor = rawSigCR ./ (1 - repmat(reshape(transpose(deadTime), [], 1, 1), 1, size(rawSigCR, 2), size(rawSigCR, 3)) * 1e-9 .* rawSigCR);
rawSigPCCor = rawSigCRCor / pc2cr;
range = (1:size(rawSigPCCor, 2)) * lData.hRes(1) + hOffset;

% �����۳�
bg = nanmean(rawSigCRCor(:, (end - 50):end, :), 2);
sigNoBg = rawSigPCCor - repmat(bg, 1, size(rawSigCRCor, 2), 1);

% �ص���������
sigNearNoBg = squeeze(sigNoBg(1, :, :));
sigFarNoBg = squeeze(sigNoBg(2, :, :));
overlapNear = interp1(ol.range, ol.overlapNear, range);
overlapFar = interp1(ol.range, ol.overlapFar, range);
sigFarNoBgCor = overlapCor(range, sigFarNoBg, overlapFar, 'glueRange', [1000, 1400]);
sigNearNoBgCor = sigNearNoBg ./ repmat(transpose(overlapNear), 1, size(sigNearNoBg, 2)) * ratioNear;

% �ź�ƴ��
sigMerge = signalMerge(sigFarNoBgCor, sigNearNoBgCor, range, [0, 100], 1, 0);

% �����
snrNear = sigNearNoBg ./ sqrt(squeeze(rawSigPCCor(1, :, :)));
snrFar = sigFarNoBg ./ sqrt(squeeze(rawSigPCCor(2, :, :)));
isNoisy = (snrNear < 1) & (snrFar < 1);

%% �״�궨
lc = getLC(lData.mTime, dbFile);
attnBsc = (sigMerge .* repmat(reshape(range, size(sigMerge, 1), 1), 1, size(sigNearNoBg, 2)).^2) / lc;
attnBsc(range <= 20, :) = NaN;

%% ����ϵ�����ݣ�ǰ���ݣ�
fExt = zeros(size(attnBsc));
feat = zeros(size(attnBsc));   % 0: no-cloud; 1: cloud; 2: unknown
lr = ones(size(attnBsc)) * 50;

[mBsc, mExt] = MolModel(range .* cos(zenithAngle / 180 * pi), 1030, 'meteor', 'standard_atmosphere');

isConverge = false;
counter = 0;
while ((~ isConverge) && (counter < 10))
    % ǰ����
    [fBsc, fExt] = quasiRetrieval(range, attnBsc, repmat(mExt', 1, size(lr, 2)), repmat(mBsc', 1, size(lr, 2)), lr, 'nIters', 10, 'flagAutoConverge', true);

    % ����ʶ��
    featBefore = feat;
    isCloud = (fExt >= cloudMinExt);
    isAerosol = (fExt < cloudMinExt);
    feat(isCloud) = 1;
    feat(isAerosol) = 0;
    feat(isNoisy) = 2;
    lr(isCloud) = 20;
    lr(isAerosol) = 45;
    lr(isNoisy) = 45;

    isSameFeatFrac = sum(sum(feat == featBefore, 2)) / numel(feat);
    isConverge = (isSameFeatFrac > 0.99);
    counter = counter + 1;
end

%% �ܼ��ȷ���
vis = ext2vis((fExt + repmat(mExt', 1, size(lr, 2))) * (1030 / 532) .^ 1.7);

% ���ݿ��ӻ�
figure('Position', [0, 20, 500, 600], 'Units', 'Pixels', 'Color', 'w');

subfig = subfigPos([0.1, 0.1, 0.83, 0.86], 3, 1, 0, 0.05);

% ���������ź�
subplot('Position', subfig(1, :), 'Units', 'normalized');
hold on;
p1 = pcolor(lData.mTime, range * 1e-3, attnBsc);
p1.EdgeColor = 'none';

xlabel('');
ylabel('���� (ǧ��)');
title('���������ź�');

xlim(tRange);
ylim([0, 10]);
caxis([0, 2e-6]);
colormap('jet');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'layer', 'top', 'Tickdir', 'out', 'XTickLabel', '', 'XTick', linspace(tRange(1), tRange(2), 5));

colorbar();

% ����ϵ��
subplot('Position', subfig(2, :), 'Units', 'normalized');
hold on;
p1 = pcolor(lData.mTime, range * 1e-3, fExt);
p1.EdgeColor = 'none';

xlabel('');
ylabel('���� (ǧ��)');
title('����ϵ�� (m-1)');

xlim(tRange);
ylim([0, 10]);
caxis([0, 1e-4]);
colormap('jet');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'layer', 'top', 'Tickdir', 'out', 'XTickLabel', '', 'XTick', linspace(tRange(1), tRange(2), 5));

colorbar();

% �ܼ���
subplot('Position', subfig(3, :), 'Units', 'normalized');
hold on;
p1 = pcolor(lData.mTime, range * 1e-3, vis);
p1.EdgeColor = 'none';

xlabel('ʱ��');
ylabel('���� (ǧ��)');
title('�ܼ��� (��)');

xlim(tRange);
ylim([0, 10]);
caxis([0, 5e4]);
colormap(gca, 'hot');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'layer', 'top', 'Tickdir', 'out', 'XTick', linspace(tRange(1), tRange(2), 5));
datetick(gca, 'x', 'HH:MM', 'Keeplimits', 'keepticks');

colorbar();

%% �ܼ��ȶԱ�

% ��ȡǰ��ɢ���ܼ���������
visDataFile = fullfile(LISMO_VARS.projectDir, 'data', 'vis35.mat');
a = load(visDataFile, 'vis35');
vis35 = a.vis35;
isInTRange = (vis35.mTime >= (tRange(1) + datenum(0, 1, 0, 8, 0, 0))) & (vis35.mTime <= (tRange(2) + datenum(0, 1, 0, 8, 0, 0)));
vis35Time = vis35.mTime(isInTRange);
vis1min = vis35.vis1min(isInTRange);
vis1minInterp = interp1(vis35Time, vis1min, lData.mTime + datenum(0, 1, 0, 8, 0, 0));
vis1minRelBias = abs(vis1minInterp - vis(idxVisCmp, :)) ./ vis1minInterp;
meanVisBias = nanmean(vis1minRelBias);

figure('Position', [0, 20, 600, 400], 'Units', 'Pixels', 'Color', 'w');

subfig = subfigPos([0.1, 0.1, 0.83, 0.8], 2, 1, 0, 0.1);

% �ܼ���ʱ������
subplot('Position', subfig(1, :), 'Units', 'normalized');

hold on;

p1 = plot(vis35Time, vis1min, 'Marker', '.', 'LineStyle', 'none', 'Color', 'k', 'DisplayName', 'ǰ��ɢ���ܼ�����');
p2 = plot(lData.mTime + datenum(0, 1, 0, 8, 0, 0), vis(idxVisCmp, :), 'Marker', 's', 'LineStyle', 'none', 'Color', 'b', 'DisplayName', '��������״�');

xlabel('');
ylabel('�ܼ��� (��)');
title('�ܼ������жԱ�');

xlim(tRange + datenum(0, 1, 0, 8, 0, 0));
ylim([0, 50000]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');
datetick(gca, 'x', 'HH:MM', 'Keepticks', 'Keeplimits');
legend([p1, p2], 'Location', 'SouthEast');


% �ܼ���ƫ��
subplot('Position', subfig(2, :), 'Units', 'normalized');

hold on;
p1 = plot(lData.mTime + datenum(0, 1, 0, 8, 0, 0), vis1minRelBias * 100, 'Color', 'b', 'DisplayName', '��������״�');
p2 = plot(lData.mTime + datenum(0, 1, 0, 8, 0, 0), ones(size(vis1minRelBias)) * 30, '--k');
p3 = plot(lData.mTime + datenum(0, 1, 0, 8, 0, 0), ones(size(vis1minRelBias)) * meanVisBias * 100, '--b', 'LineWidth', 2);
hold off;

text(0.6, 0.8, sprintf('ƽ�����ƫ��Ϊ%5.2f%%', meanVisBias * 100), 'Units', 'Normalized');

xlabel('');
ylabel('���ƫ�� (%)');

xlim(tRange + datenum(0, 1, 0, 8, 0, 0));
ylim([0, 50]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');
datetick(gca, 'x', 'HH:MM', 'Keepticks', 'Keeplimits');