% 批处理能见度激光雷达数据，反演得到消光系数和能见度
%
% 作者：殷振平
% 邮箱: zp.yin@whu.edu.cn
% 日期：2026-05-29

%% Parameter Definition
l0Folder = 'G:\backup\vis-lidar\20260428-tianjin-evaluation\前散仪对比数据\VIS1\L0';   % 雷达原始数据文件目录
l1Folder = 'G:\backup\vis-lidar\20260428-tianjin-evaluation\前散仪对比数据\VIS1\L1';   % 雷达产品文件目录（用于输出结果存储）
savePath = 'H:\research\vislidar-intercomparison\tianjing\VIS1';   % 输出结果目录
overlapFile = 'VIS1_overlap_factor.txt';   % 重叠因子文件，重叠因子通过水平方法计算
lr = 50;   % 雷达比
refH = [8, 9];   % Fernald方法参考高度（千米）
AEConvFactor = (1030/550) ^1;   % 波长指数（转换因子）
flagOLCor = false;   % 是否进行重叠因子修正
flagSaveProfileFiles = false;   % 是否保存单个剖面数据文件
eleAgl = 6;   % 雷达仰角（度）
flagDisplay = 'on';

%% Load Data

% read overlap factor
olHeight = [];
olVal = [];
if exist(fullfile(savePath, overlapFile), 'file') == 2
    fid = fopen(overlapFile, 'r');

    tmp = textscan(fid, '%f%f', 'HeaderLines', 1, 'Delimiter', ' ');
    olHeight = tmp{1};
    olVal = tmp{2};

    fclose(fid);
end

%% find dates
dateFolders = listdir(l0Folder, '\w*', 1);

for iFolder = 9:length(dateFolders)

    fprintf('Finished %6.2f%%: processing %s\n', (iFolder - 1) / length(dateFolders) * 100, dateFolders{iFolder});

    close all;

    thisDate = datenum(basename(dateFolders{iFolder}), 'yyyymmdd');

    %% 加载原始数据
    %（已经从dat文件转换成了mat格式文件），转换脚本为convert_lidar_data_2_mat.m
    l0DataFile = fullfile(l0Folder, sprintf('%s_vis_lidar_l0.mat', datestr(thisDate, 'yyyy-mm-dd')));
    load(l0DataFile);

     %% Rayleigh scattering
    [mBsc, mExt] = MolModel(range * sin(eleAgl / 180 * pi), 1030, 'meteor', 'standard_atmosphere');
    [~, mExt550] = MolModel(range * sin(eleAgl / 180 * pi), 550, 'meteor', 'standard_atmosphere');

    %% 数据预处理
    thisBG = nanmean(lidarSig((end - 70):(end - 10), :), 1);
    lidarSigCor = lidarSig - repmat(thisBG, size(thisBG, 1), 1);
    if flagOLCor
        ol = interp1(olHeight, olVal, range);
        lidarSigCor = lidarSigCor ./ repmat(ol, 1, size(lidarSigCor, 2));
    end
    lidarRCS = lidarSigCor .* repmat(range, 1, size(lidarSigCor, 2)).^2;

    %% 可视化

    % 距离修正信号时空高度图
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w', 'visible', flagDisplay);

    hold on;
    p1 = pcolor(mTime, range * 1e-3, lidarRCS);
    p1.EdgeColor = 'none';
    hold off;

    xlabel('时间');
    ylabel('距离 (千米)');
    title('距离修正信号');

    xlim([min(mTime), max(mTime)]);
    ylim([0, 10]);
    caxis([0, 5e10]);

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'layer', 'Top', 'TickDir', 'out', 'FontSize', 11);

    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    colorbar();
    export_fig(gcf, fullfile(savePath, sprintf('%s_RCS.png', datestr(mTime(1), 'yyyy-mm-dd'))), '-r300');

    %% extinction
    extMat_Xian = NaN(size(lidarSigCor));
    extMat_Fernald = NaN(size(lidarSigCor));
    for iPrf = 1:size(lidarSigCor, 2)
        ext1 = extRet_Fernald(range, lidarSigCor(:, iPrf), thisBG(iPrf), mBsc, mExt, 'minSNR', 3, 'hFullOL', 330);
        ext3 = extRet_Xian(range, lidarSigCor(:, iPrf), thisBG(iPrf), 'minSNR', 3, 'rangeFullOverlap', 330);

        extMat_Xian(:, iPrf) = ext3;
        extMat_Fernald(:, iPrf) = ext1;
    end

    %% extinction to visibility
    visMat_Xian = ext2vis(extMat_Xian);
    visMat_Fernald = ext2vis(extMat_Fernald * AEConvFactor + mExt550);

    %% Save to .mat file
    save(fullfile(l1Folder, sprintf('%s_vis_lidar_l1_exp.mat', datestr(mTime(1), 'yyyy-mm-dd'))), ...
        'mTime', 'range', 'visMat_Xian', 'extMat_Xian', 'visMat_Fernald', 'extMat_Fernald');

    %% Save to ASCII file
    if flagSaveProfileFiles

        subSaveFolder = fullfile(savePath, datestr(thisDate, 'yyyymmdd'));
        if exist(subSaveFolder, 'dir') ~= 7
            mkdir(subSaveFolder);
        end

        rangeQC = range(range >= 0);

        for iPrf = 1:size(extMat_Xian, 2)

            % quality control
            extMatQC = extMat_Xian(range > 0, iPrf) * 1e3;
            extMatQC(extMatQC <= 0) = -9999;
            extMatQC(isnan(extMatQC)) = -9999;
            visMatQC = visMat_Xian(range > 0, iPrf) * 1e-3;
            visMatQC(visMatQC <= 0) = -9999;
            visMatQC(visMatQC > 50) = 50;
            visMatQC(isnan(visMatQC)) = -9999;

            % save files
            outputFile = fullfile(subSaveFolder, sprintf('vis_lidar_%s_level1.txt', datestr(mTime(iPrf), 'yyyymmdd_HHMMSS')));

            fid = fopen(outputFile, 'w');
            fprintf(fid, 'instrument: visibility lidar-1030nm\n');
            fprintf(fid, 'measurement time: %s\n', datestr(mTime(iPrf), 'yyyy-mm-dd HH:MM:SS'));
            fprintf(fid, 'elevation angle (degree): %.1f\n', eleAgl);
            fprintf(fid, 'azimuth angle (degree): 0.0\n');
            fprintf(fid, 'geolocation (lat, lon): 30.5, 114.5\n');
            fprintf(fid, '\n');
            fprintf(fid, 'Height(m); 消光系数(km-1); 能见度(km)\n');
            for iH = 1:length(rangeQC)
                fprintf(fid, '%.1f; %.5f; %.2f\n', rangeQC(iH), extMatQC(iH), visMatQC(iH));
            end

            fclose(fid);
        end
    end

    %% Display

    %% Xian
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w', 'visible', flagDisplay);

    subplot('Position', [0.1, 0.2, 0.72, 0.7], 'Units', 'Normalized');

    hold on;
    p1 = pcolor(mTime, range * 1e-3, visMat_Xian * 1e-3);
    p1.EdgeColor = 'none';

    xlabel('时间');
    ylabel('距离 (千米)');
    title(sprintf('%s 能见度 (大舜算法)', datestr(mTime(1), 'yyyy-mm-dd')));

    xlim([min(mTime), max(mTime)]);
    ylim([0, 10]);
    caxis([0, 20]);

    colormap('jet');

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'XTick', linspace(min(mTime), max(mTime), 5), 'layer', 'Top', 'TickDir', 'out', 'FontSize', 11);

    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    cb = colorbar('Position', [0.85, 0.3, 0.03, 0.5], 'Units', 'Normalized');
    titleHandle = get(cb, 'Title');
    set(titleHandle, 'String', '[千米]');
    set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

    export_fig(gcf, fullfile(savePath, sprintf('%s_vis_thi.png', datestr(mTime(1), 'yyyy-mm-dd'))), '-r300');

    %% Fernald
    figure('Position', [0, 30, 450, 250], 'Units', 'Pixels', 'Color', 'w', 'visible', flagDisplay);

    subplot('Position', [0.1, 0.2, 0.72, 0.7], 'Units', 'Normalized');

    hold on;
    p1 = pcolor(mTime, range * 1e-3, visMat_Fernald * 1e-3);
    p1.EdgeColor = 'none';

    xlabel('时间');
    ylabel('距离 (千米)');
    title(sprintf('%s 能见度 (Fernald算法)', datestr(mTime(1), 'yyyy-mm-dd')));

    xlim([min(mTime), max(mTime)]);
    ylim([0, 10]);
    caxis([0, 20]);

    colormap('jet');

    set(gca, 'YMinorTick', 'on', 'Box', 'on', 'XTick', linspace(min(mTime), max(mTime), 5), 'layer', 'Top', 'TickDir', 'out', 'FontSize', 11);

    datetick(gca, 'x', 'HH:MM', 'keeplimits', 'keepticks');

    cb = colorbar('Position', [0.85, 0.3, 0.03, 0.5], 'Units', 'Normalized');
    titleHandle = get(cb, 'Title');
    set(titleHandle, 'String', '[千米]');
    set(cb, 'TickDir', 'in', 'Box', 'on', 'TickLength', 0.02);

    export_fig(gcf, fullfile(savePath, sprintf('%s_vis_thi_Fernald.png', datestr(mTime(1), 'yyyy-mm-dd'))), '-r300');

end