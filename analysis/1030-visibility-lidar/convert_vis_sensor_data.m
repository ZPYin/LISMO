% convert visibility sensor data.
% 
% Author: Zhenping Yin
% Date: 2026-05-27
% Email: zp.yin@whu.edu.cn

%% Parameter Definition
% visFile = 'G:\backup\vis-lidar\20260428-tianjin-evaluation\前散仪对比数据\2026年05月01日00时+至+2026年05月07日00时+观测试验基地+____+能见度对比+分钟数据.xls';
visFile = 'G:\backup\vis-lidar\20251221-wuhan\25-12-21.log';

%% Read Data
[~, ~, ext] = fileparts(visFile);

visData = struct();
visData.mTime = [];
visData.vis = [];
switch ext
    case {'.xls', '.xlsx'}
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
        fid = fopen(visFile, 'r');

        lineCounter = 0;
        while ~feof(fid)
            fprintf('Line %d\n', lineCounter);

            thisLine1 = fgetl(fid);
            lineCounter = lineCounter + 1;
            if (~contains(thisLine1, 'RECV ASCII'))
                continue;
            else
                thisLine2 = fgetl(fid);
                lineCounter = lineCounter + 1;
                if (~contains(thisLine2, 'PW')) || (~contains(thisLine2, '/////'))
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

% legend(p1, 'Location', 'NorthWest');

export_fig(gcf, 'vis_sensor_overview.png', '-r300');