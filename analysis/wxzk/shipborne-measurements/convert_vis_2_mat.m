% 将前向散射能见度仪测量的能见度时间序列文件转换为MAT数据

clc;
close all;

visFile = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\海试数据、效果图\Visibility.txt';
% visFile = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\岸基前向散射能见度仪\Visibility-岸基.txt';
saveFile = 'visibility.mat';

%% 读取数据文件
mTime = [];
vis = [];
fid = fopen(visFile, 'r');

while ~feof(fid)
    thisLine = fgetl(fid);

    strs = strsplit(thisLine, ',');

    mTime = cat(2, mTime, datenum(strs{1}, 'yyyy-mm-dd HH:MM:SS'));
    vis = cat(2, vis, str2double(strs{2}));
end

fclose(fid);

save(saveFile, 'mTime', 'vis');