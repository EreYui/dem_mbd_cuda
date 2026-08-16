%% 根据刚体位姿输出结果，生成一系列刚体vtk，用于paraview可视化

clear
clc
parameter;
%% 参数
bodyPath = [Project_filefolder,'Data\',datafilefolder,'\InputFile\BodySet\'];
bodyvtkPath=[Project_filefolder,'Data\',datafilefolder,'\OutputFile\vtk_bodys\'];
motionDataPath = [Project_filefolder,'Data\',datafilefolder,'\OutputFile\state_bodys\'];

original_vtk_paths={
    [bodyPath,'cube_150mm_surface.vtk']%,
    % [bodyPath,'plate2.vtk'],
    %[bodyPath,'plate3.vtk']
    };

output_dir = {
    [bodyvtkPath,'body1']%,
    % [bodyvtkPath,'body2'],
    %[bodyvtkPath,'body3']
    };

nfile = file_num_body; % 文件数量
nbody = num_body;      % 刚体数量
count = StepOutput;    % 输出文件步长间隔
%% 文件转换
pos = zeros(nbody, nfile, 3);
quat = zeros(nbody, nfile, 4);
for timeStep = 1:nfile%
    % 生成文件名
    resfile = fullfile(motionDataPath, ['At.', sprintf('%08d', timeStep * count), '.bt']);
    tmp = load(resfile);
    
    % 处理每个刚体数据
    for bodyID = 1:nbody
        row_num = bodyID;%(bodyID-1)*count+1;  % 计算目标行号
        line_data = tmp(row_num,:); % 解析行数据
        
        % 提取所需数据
        pos(bodyID,timeStep,:) = [line_data(2), line_data(3), line_data(4)];
        quat(bodyID,timeStep,:) =[line_data(11), line_data(12), line_data(13), line_data(14)];
    end
    fprintf('正在处理: %d/%d\n', timeStep, nfile);
end
generateDynamicVTKs(original_vtk_paths, pos, quat, output_dir);
disp('刚体（系）vtk已输出。');

function generateDynamicVTKs(original_vtk_paths, pos, quat, output_dir)
% 生成各时刻下各刚体的VTK文件
% 输入参数：
%   original_vtk_paths: 单元格数组，每个元素是刚体的原始VTK文件路径
%   pos: N×n×3 数组，pos(i,t,:)表示刚体i在时刻t的位置 [x,y,z]
%   quat: N×n×4 数组，quat(i,t,:)表示刚体i在时刻t的四元数，格式为[w,x,y,z]
%   output_dir: 输出目录

N = numel(original_vtk_paths); % 刚体数量
n = size(pos, 2); % 时间步数

% 预加载所有刚体的VTK数据
headers = cell(N, 1);
points_data = cell(N, 1);
remaining_lines = cell(N, 1);
for i = 1:N
    [headers{i}, points_data{i}, remaining_lines{i}] = read_original_vtk(original_vtk_paths{i});
    points_data{i} = points_data{i};
end

% 确保输出目录存在
for i=1:N
    if ~exist(output_dir{i}, 'dir')
        mkdir(output_dir{i});
    end
end

% 遍历每个时间步和刚体
for t = 1:n
    for i = 1:N
        % 当前位姿
        current_pos = squeeze(pos(i,t,:))';
        current_quat = squeeze(quat(i,t,:))';

        % 计算旋转矩阵（四元数格式需为[w,x,y,z]）
        quat_normalized = current_quat / norm(current_quat);
        R = quat2rotm(quat_normalized);

        % 变换节点坐标
        transformed_points = (R * points_data{i}')' + current_pos;

        % 生成新的节点坐标行
        num_points = size(transformed_points, 1);
        new_point_lines = cell(num_points, 1);
        for j = 1:num_points
            new_point_lines{j} = sprintf('%.9e %.9e %.9e', ...
                transformed_points(j,1), transformed_points(j,2), transformed_points(j,3));
        end

        % 构建新VTK内容
        new_vtk_content = [headers{i}; new_point_lines; remaining_lines{i}];

        % 写入文件
        output_filename = fullfile(output_dir{i}, ...
            sprintf('body%d_time%08d.vtk', i, t));
        fid = fopen(output_filename, 'w', 'n', 'US-ASCII');
        fprintf(fid, '%s\n', new_vtk_content{:});
        fclose(fid);
    end
end
end

function [header, points_data, remaining_lines] = read_original_vtk(filename)
% 读取原始VTK文件，分离头部、节点数据和剩余内容
lines = readlines(filename);
point_line_idx = find(startsWith(lines, 'POINTS'), 1);
if isempty(point_line_idx)
    error('POINTS line not found in VTK file');
end

% 解析节点数量
parts = strsplit(lines{point_line_idx});
num_points = str2double(parts{2});

% 提取头部（包含POINTS行）
header = lines(1:point_line_idx);

% 提取节点数据
points_data = zeros(num_points, 3);
for j = 1:num_points
    points_data(j,:) = sscanf(lines{point_line_idx + j}, '%f %f %f')';
end

points_data = points_data*0.001;
% 剩余行
remaining_lines = lines(point_line_idx + num_points + 1 : end);
end