% filePath: 原始文件的完整路径（例如 'F:\CoupDyn\Data\DATA 5\InputFile\BodySet\leg.bt'）
% copies: 需要复制的份数（例如 12）
filePath = 'F:\CoupDyn\Data\DATA\InputFile\BodySet\leg.vtk';
copies = 12;

% 检查文件是否存在
if ~isfile(filePath)
    error('原始文件不存在: %s', filePath);
end

% 分解文件路径信息
[folder, baseName, ext] = fileparts(filePath);

% 检查文件夹是否存在，不存在则创建
if ~isfolder(folder)
    mkdir(folder);
    fprintf('已创建目录: %s\n', folder);
end

% 复制文件
for i = 1:copies
    newFileName = sprintf('%s%d%s', baseName, i, ext);
    newFilePath = fullfile(folder, newFileName);

    % 复制文件（若目标已存在则覆盖）
    copySuccess = copyfile(filePath, newFilePath);

    if copySuccess
        fprintf('已创建: %s\n', newFilePath);
    else
        warning('复制失败: %s', newFilePath);
    end
end

fprintf('操作完成! 共创建 %d 个副本文件\n', copies);