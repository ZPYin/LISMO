% 将前向散射能见度仪测量的能见度时间序列文件转换为MAT数据
%
% 殷振平
% 2024-01-21

clc;
close all;

%% 定义参数
visFile = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\岸基前向散射能见度仪\Visibility-岸基.txt';
saveFile = 'vis-fixed-platform.mat';

%% 读取数据文件
fid = fopen(visFile, 'r');

data = textscan(fid, '%s%f', 'delimiter', ',', 'HeaderLines', 0);

fclose(fid);

mTime = NaN(1, length(data{1}));
vis = NaN(1, length(data{1}));
for iLine = 1:length(data{1})
    mTime(iLine) = datenum(data{1}{iLine}, 'yyyy-mm-dd HH:MM:SS');
    vis(iLine) = data{2}(iLine);
end

save(saveFile, 'mTime', 'vis');