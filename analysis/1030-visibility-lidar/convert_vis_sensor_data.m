% convert visibility sensor data.
% 
% Author: Zhenping Yin
% Email: zp.yin@whu.edu.cn
% Date: 2026-05-27

parentPath = fileparts(mfilename('fullpath'));

%% Parameter Definition
% visFile = 'G:\backup\vis-lidar\20260428-tianjin-evaluation\前散仪对比数据\2026年05月01日00时+至+2026年05月07日00时+观测试验基地+____+能见度对比+分钟数据.xls';
visFile = 'G:\backup\vis-lidar\20251221-wuhan\25-12-21.log';  % 目前支持两种格式.xls/.xlsx和.log，其他格式可以参考这两种进行修改

%% Read Data
[~, ~, ext] = fileparts(visFile);

visData = struct();
visData.mTime = [];
visData.vis = [];

switch ext
    case {'.xls', '.xlsx'}
        % Excel格式文件，这种文件主要是应用于天津标校场的前向散射能见度仪输出结果
        [~, ~, dataCell] = xlsread(visFile);

        %% Data Conversion
        nRec = size(dataCell, 1) - 1;
        mTime = NaN(1, nRec);
        vis = NaN(1, nRec);
        for iR = 1:nRec
            mTime(iR) = datenum(dataCell{iR + 1, 2}, 'yyyy/mm/dd HH:MM:SS');
            vis(iR) = dataCell{iR + 1, 4};
        end
        visData.mTime = mTime;
        visData.vis = vis;

    case '.log'
        % 文本格式文件，主要是武汉量子院前向散射能见度仪输出文件
        % 其实正常前向散射能见度仪应该输出.csv格式，但是有人把log文件拷过来了
        % 硬着头皮处理而已
        fid = fopen(visFile, 'r');

        lineCounter = 0;
        while ~feof(fid)
            fprintf('Line %d\n', lineCounter);

            thisLine1 = fgetl(fid);
            lineCounter = lineCounter + 1;

            if (~contains(thisLine1, 'RECV ASCII'))
                % 这一部分主要是做一些数据清洗，去掉一些无效行
                continue;
            else
                thisLine2 = fgetl(fid);
                lineCounter = lineCounter + 1;
                if (~contains(thisLine2, 'PW')) || (~contains(thisLine2, '/////'))
                    % 这一部分主要是做一些数据清洗，去掉一些无效行
                    continue;
                else
                    thisMTime = datenum(thisLine1(2:24), 'yyyy-mm-dd HH:MM:SS.FFF');
                    thisVis = str2double(thisLine2(10:15));

                    visData.mTime = cat(2, visData.mTime, thisMTime);
                    visData.vis = cat(2, visData.vis, thisVis);
                end
            end
        end

        fclose(fid);

    otherwise
        error('Unsupported file format: %s', ext);
end

%% Save Data
save('vis-sensor-data.mat', 'visData');

%% Display
figure('Position', [0, 30, 500, 250], 'Units', 'Pixels', 'Color', 'w');

hold on;
p1 = plot(visData.mTime, visData.vis * 1e-3, '-r');
hold off;

xlim([min(visData.mTime), max(visData.mTime)]);
ylim([0, 50]);

xlabel('时间');
ylabel('能见度 (千米)');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 11);
datetick(gca, 'x', 'mm-dd', 'keeplimits', 'keepticks');

export_fig(gcf, fullfile(parentPath, 'vis_sensor_overview.png'), '-r300');