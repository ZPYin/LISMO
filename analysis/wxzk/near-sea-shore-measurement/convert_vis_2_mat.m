% ��ǰ��ɢ���ܼ����ǲ������ܼ���ʱ�������ļ�ת��ΪMAT����
%
% ����ƽ
% 2024-01-21

clc;
close all;

%% �������
visFile = 'C:\Users\ZPYin\Documents\Data\wxzk_fog_measurements\����ǰ��ɢ���ܼ�����\Visibility-����.txt';
saveFile = 'vis-fixed-platform.mat';

%% ��ȡ�����ļ�
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