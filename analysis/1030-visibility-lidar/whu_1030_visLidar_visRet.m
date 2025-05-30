% ���ݵõ��人��ѧ�ܼ��ȼ����״��Ʒ
% - ���ܽ�����ɢ��ϵ��
% - ���ܽ�����ϵ��
% - �ܼ���
% ���ߣ�����ƽ
% ���䣺zp.yin@whu.edu.cn
% ���ڣ�2025��5��21��

clc;
close all;

%% Parameter Definition
dataFile = 'C:\Users\zhenp\OneDrive\Desktop\VisibilityInversion (3)\VisibilityInversion (3)\VisibilityInversion\S001-NAP001-Test CMA-017-210927-110120.dat';   % �ź�����
resPath = 'C:\Users\zhenp\OneDrive\Desktop';   % ��Ʒ���Ŀ¼
figPath = 'C:\Users\zhenp\OneDrive\Desktop';   % ͼƬ������Ŀ¼
flagOLCor = false;   % �Ƿ�����ص���������
olFile = 'overlap_20250101.mat';   % �ص������ļ�
AEConvFactor = (1030/550) ^ 1;   % ���ܽ�����ϵ������ת������
hFullOL = 500;   % ��ȫ�����ӳ��߶ȣ�m��
lidarRatio = 50;   % �����״�� (sr)
distOffset = -55;   % Ԥ�������������ͨ���źŵ�һ����ֵ�����жϣ�
visible = 'on';   % �Ƿ���н�����ӻ�

%% Read Data
if exist(dataFile, 'file')
    thisData = readALADat(dataFile);
else
    warning('Data file does not exist: %s', dataFile);
    return;
end

%% Preprocessing
range = transpose(((1:thisData.nBins) - 0.5 + distOffset) * thisData.hRes(1));
bg = mean(thisData.rawSignal((end - 30):(end - 5), :), 1);
noise = std(thisData.rawSignal((end - 30):(end - 5), :), 0, 1);
signal = thisData.rawSignal - repmat(bg, thisData.nBins, 1);
angle = thisData.ele;
mTime = thisData.mTime;
if flagOLCor
    ol = load(olFile);
    signal = signal ./ repmat(transpose(ol.ov), 1, size(thisData.rawSignal, 2));
end
rcs = signal .* repmat(range, 1, size(thisData.rawSignal, 2)).^2;   % ���������ź�
snr = signal ./ repmat(noise, thisData.nBins, 1);   % �����
height = range .* sin(angle / 180 * pi);   % �߶ȣ�m��

%% Rayleigh Scattering
[mBsc, mExt] = MolModel(height, 1030, 'meteor', 'standard_atmosphere');
[~, mExt550] = MolModel(height, 550, 'meteor', 'standard_atmosphere');

%% extinction retrieval
extMat_Fernald = extRet_Fernald(range, signal, bg, mBsc, mExt, ...
    'snr', snr, ...
    'minSNR', 3, ...
    'hFullOL', 500, ...
    'lr', lidarRatio);
bscMat_Fernald = extMat_Fernald ./ lidarRatio;

%% Extinction to Visibility
visMat_Fernald = ext2vis(extMat_Fernald * AEConvFactor + mExt550);

%% Save Results
fid = fopen(fullfile(resPath, 'test_output.txt'), 'w');

fprintf(fid, '������Ϣ: \n');
fprintf(fid, '�ļ���: %s\n', basename(dataFile));

fprintf(fid, '\n���ݽ��\n');
fprintf(fid, '����(m); ���������ź�; �����; ���ܽ�����ɢ��ϵ�� (km-1sr-1); ���ܽ�����ϵ�� (km-1); �ܼ��� (km);\n');

for iBin = 1:thisData.nBins
    fprintf(fid, '%f; %f; %f; %f; %f; %f;\n', range(iBin), rcs(iBin), snr(iBin), bscMat_Fernald(iBin) * 1e3, extMat_Fernald(iBin) * 1e3, visMat_Fernald(iBin) * 1e-3);
end

fclose(fid);

%% Display

% range corrected signal
figure('Position', [100, 100, 600, 600], 'color', 'w', 'visible', visible);

% range corrected signal
subplot(311);

hold on;
YRcs = rcs;
YRcs(YRcs <= 0) = NaN;
plot(range * 1e-3, YRcs, 'b');
hold off;

xlabel('���� (km)');
ylabel('���������ź� (a.u.)');

xlim([0, 15]);
ylim('auto');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'YTick', 10.^(6:14), ...
    'Box', 'on', 'FontSize', 12, 'YScale', 'log');

% aerosol backscatter and extinction
subplot(312);

hold on;
p1 = plot(range * 1e-3, extMat_Fernald * 1e3, 'b', 'DisplayName', '���ܽ�');
p2 = plot(range * 1e-3, mExt * 1e3, 'r', 'DisplayName', '��������');
hold off;

xlabel('���� (km)');
ylabel('����ϵ�� (km-1)');

xlim([0, 15]);
ylim([-0.005, 0.1]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'YTick', 0:0.2:1, ...
    'Box', 'on', 'FontSize', 12, 'YScale', 'Linear');

l = legend([p1, p2], 'Location', 'NorthEast');
l.FontSize = 11;

% visibility
subplot(313);

hold on;
plot(range * 1e-3, visMat_Fernald * 1e-3, 'b', 'DisplayName', '���ܽ�');
plot([-100, 100], [50, 50], '-.k');
hold off;

xlabel('���� (km)');
ylabel('�ܼ��� (km)');

xlim([0, 15]);
ylim([0, 60]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'YTick', 0:10:50, ...
    'Box', 'on', 'FontSize', 12, 'YScale', 'Linear');

export_fig(gcf, fullfile(figPath, 'test_results_output.png'), '-r300');
