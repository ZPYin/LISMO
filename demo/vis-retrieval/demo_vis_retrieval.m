% 用于测试能见度反演算法的样例
%
% 作者：殷振平
% 邮箱：zp.yin@whu.edu.cn
% 日期：2026-05-29

clc;
close all;
parentPath = fileparts(mfilename('fullpath'));

%% 参数定义
lidarDataFile = fullfile(parentPath, 'inputs', 'Z_RADA_I_51544_20260503000030_O_VISLIDAR_EP-VL10_L0.dat');
olFile = fullfile(parentPath, 'inputs', 'VIS_xglz_overlap_factor.txt');
nPretriggers = 30;

%% 数据读取

% 读取重叠因子
olHeight = [];
olVal = [];
if exist(olFile, 'file') == 2
    fid = fopen(olFile, 'r');

    tmp = textscan(fid, '%f%f', 'HeaderLines', 1, 'Delimiter', ' ');
    olHeight = tmp{1};
    olVal = tmp{2};

    fclose(fid);
end

% 读取激光雷达数据 (1 channel)
lidarData = readALADat(lidarDataFile, 'datatype', 2);

%% 信号预处理

range = ((1:lidarData.nBins) - nPretriggers + 0.5) * lidarData.hRes(1);

bg = nanmean(lidarData.rawSignal((end - 40):(end - 10)));   % 背景噪声水平
signal = lidarData.rawSignal - bg;   % 去除背景噪声

ol = ones(size(lidarData.rawSignal));
if ~ isempty(olHeight)
    ol = interp1(olHeight, olVal, range);
end
sigCor = transpose(signal) ./ ol;   % 重叠校正

%% 能见度反演

% 消光系数
tExt = extRet_Xian(range, sigCor, bg, 'rangeFullOverlap', 500.0, 'minSNR', 2);

% 后向散射系数
tBsc = tExt / 50;

% 能见度
vis = ext2vis(tExt);

%% 保存结果
rangeQC = range(range >= 0);

% quality control
extQC = tExt(range > 0) * 1e3;
extQC(extQC <= 0) = -9999;
extQC(isnan(extQC)) = -9999;
visQC = vis(range > 0) * 1e-3;
visQC(visQC <= 0) = -9999;
visQC(visQC > 50) = 50;
visQC(isnan(visQC)) = -9999;

% save files
outputFile = fullfile(parentPath, 'outputs', sprintf('vis_lidar_%s_level1.txt', datestr(lidarData.mTime, 'yyyymmdd_HHMMSS')));

fid = fopen(outputFile, 'w');
fprintf(fid, 'instrument: visibility lidar-1030nm\n');
fprintf(fid, 'measurement time: %s\n', datestr(lidarData.mTime, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '\n');
fprintf(fid, 'Height(m); 消光系数(km-1); 能见度(km)\n');
for iH = 1:length(rangeQC)
    fprintf(fid, '%.1f; %.5f; %.2f\n', rangeQC(iH), extQC(iH), visQC(iH));
end

fclose(fid);