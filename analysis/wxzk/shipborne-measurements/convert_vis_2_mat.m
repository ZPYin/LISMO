% ��ǰ��ɢ���ܼ����ǲ������ܼ���ʱ�������ļ�ת��ΪMAT����

clc;
close all;

visFile = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\�������ݡ�Ч��ͼ\Visibility.txt';
% visFile = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\����ǰ��ɢ���ܼ�����\Visibility-����.txt';
saveFile = 'visibility.mat';

%% ��ȡ�����ļ�
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