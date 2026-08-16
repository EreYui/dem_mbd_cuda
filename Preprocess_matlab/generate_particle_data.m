%generate_particle_data()
%% 从已有数据中按需求进行提取，编号重新编，从0开始逐渐增加
% 1.按空间位置，将特定位置的提取出来
% 2.将颗粒位置进行平移
% 3.将颗粒块进行复制平移，组成更大范围内的颗粒系
% 4.对颗粒进行缩放，相应地质量和空间位置也会有变化
%% 数据处理
% 参数设置
p.mode = 3;        % 工作模式：1-提取 2-平移 3-复制组合
p.dxyz = [0.5 0 0];  % 平移量
p.xyz = [-0.5,0.5,-0.5,0.5,0,0.5]; % 目标范围[xmin,xmax,ymin,ymax,zmin,zmax]
p.z = 0;           % 是否考虑z轴范围
p.r = 1;           % 缩放倍数
%% 文件路径
ptfile = 'F:\CoupDyn\Preprocess_matlab\已有的颗粒文件\pt3.bt';
xyz = [0,0.5,0,0.5,0,0.5]; % 原始数据范围
% 读取数据
data = load(ptfile);
%% 数据缩放处理
if p.r ~= 1
    data(:,3) = data(:,3) * (p.r^3);   % 质量缩放
    data(:,4) = data(:,4) * p.r;       % 半径缩放
    data(:,5:7) = data(:,5:7) * p.r;   % 位置缩放
    data(:,8:10) = data(:,8:10) * p.r; % 速度缩放
end
%% 主处理逻辑
switch p.mode
    case 1  % 模式1：范围提取
        new_data = extract_particles(data, p.xyz, p.z);
        
    case 2  % 模式2：数据平移
        data = translate_particles(data, p.dxyz);
        new_data = data;
        
    case 3  % 模式3：复制组合
        % 计算原始块尺寸
        dx = xyz(2)-xyz(1);
        dy = xyz(4)-xyz(3);
        dz = xyz(6)-xyz(5);
        
        % 计算复制次数
        nx = ceil((p.xyz(2)-p.xyz(1))/dx);
        ny = ceil((p.xyz(4)-p.xyz(3))/dy);
        nz = p.z*ceil((p.xyz(6)-p.xyz(5))/dz) + ~p.z;
        
        % 将原数据平移到目标位置的起始点
        trans = [p.xyz(1)-xyz(1) p.xyz(3)-xyz(3) p.xyz(5)-xyz(5) ];
        if norm(trans) ~=0
            data = translate_particles(data, trans);
        end
        
        % 生成所有平移组合
        [X,Y,Z] = meshgrid(0:nx-1, 0:ny-1, 0:nz-1);
        translations = [X(:)*dx, Y(:)*dy, Z(:)*dz];
        
        % 复制和平移数据
        new_data = [];
        for i = 1:size(translations,1)
            translated = translate_particles(data, translations(i,:));
            new_data = [new_data; translated];
        end
        
        % 范围提取
        new_data = extract_particles(new_data, p.xyz, p.z);
end

% 重新编号并保存
new_data(:,1) = 0:size(new_data,1)-1;
save_particles(new_data, 'output.txt');


% 范围提取函数
function extracted = extract_particles(data, range, zflag)
xmin = range(1); xmax = range(2);
ymin = range(3); ymax = range(4);
zmin = range(5); zmax = range(6);

mask = (data(:,5)-data(:,4)) >= xmin & (data(:,5)+data(:,4)) <= xmax & ...
    (data(:,6)- data(:,4)) >= ymin & (data(:,6) + data(:,4))<= ymax;

if zflag
    mask = mask & (data(:,7)-data(:,4)) >= zmin & (data(:,7)+data(:,4)) <= zmax;
end
extracted = data(mask,:);
end

% 数据平移函数
function translated = translate_particles(data, dxyz)
translated = data;
translated(:,5) = translated(:,5) + dxyz(1); % x
translated(:,6) = translated(:,6) + dxyz(2); % y
translated(:,7) = translated(:,7) + dxyz(3); % z
end

% 数据保存函数
function save_particles(data, filename)
fid = fopen(filename, 'w');
for i = 1:size(data,1)
    fprintf(fid, ['%d %d %.12f %.12f '...     % 编号、状态、质量、半径
        '%.12f %.12f %.12f '...               % 位置
        '%.12f %.12f %.12f '...               % 速度
        '%.12f %.12f %.12f %.12f '...         % 四元数
        '%.12f %.12f %.12f\n'], data(i,:));   % 角速度
end
fclose(fid);
end