
parameter;
targetDir = [Project_filefolder,'Data\',datafilefolder,'\OutputFile'];

clearFolderContent(targetDir);

function clearFolderContent(targetDir)
    % 清空目标文件夹下所有文件（包括子文件夹），保留文件夹结构
    % 获取所有文件的路径（递归搜索）
    allFiles = dir(fullfile(targetDir, '**', '*.*'));
    
    % 排除文件夹条目（保留纯文件）
    allFiles = allFiles(~[allFiles.isdir]);
    
    % 检查是否找到文件
    if isempty(allFiles)
        disp('目标文件夹中没有发现任何文件。');
        return;
    end
    
    % 逐文件删除
    for i = 1:length(allFiles)
        filePath = fullfile(allFiles(i).folder, allFiles(i).name);
        try
            delete(filePath);
        catch ME
            warning('无法删除文件: %s\n原因: %s', filePath, ME.message);
        end
    end
    
    % 完成报告
    fprintf('已删除 %d 个文件，保留所有文件夹结构。\n', length(allFiles));
end