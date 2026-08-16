%% 判断外法线方向是否正确


% 读取.bt文件
filename = 'leg.bt';
fid = fopen(filename, 'r');

% 读取节点和三角形数量
n_nodes = fscanf(fid, '%d', 1);
n_tris = fscanf(fid, '%d', 1);

% 读取节点坐标（n×3矩阵）
nodes = fscanf(fid, '%f %f %f', [3, n_nodes])';
nodes = nodes(:, 1:3); % 确保格式正确

% 读取三角形索引（n×3矩阵）
tris = fscanf(fid, '%d %d %d', [3, n_tris])';

fclose(fid);

% 预分配结果存储
normals = zeros(n_tris, 3);
centroids = zeros(n_tris, 3);
is_outward = false(n_tris, 1);

% 遍历所有三角形
for i = 1:n_tris
    % 获取三个顶点索引（注意索引从1开始）
    idx = tris(i, :);
    v1 = nodes(idx(1), :);
    v2 = nodes(idx(2), :);
    v3 = nodes(idx(3), :);
    
    % 计算边向量
    a1 = v2 - v1;
    a2 = v3 - v1;
    
    % 计算法线（叉乘）
    normal = cross(a1, a2);
    normal = normal / norm(normal); % 单位化
    
    % 计算三角形中心
    centroid = mean([v1; v2; v3]);
    
    % 判断法线方向是否朝外（假设球心在原点）
    center_to_surface = centroid - [0 0 0]; % 球心到表面的方向
    dot_product = dot(normal, center_to_surface);
    
    % 存储结果
    normals(i, :) = normal;
    centroids(i, :) = centroid;
    is_outward(i) = dot_product > 0;
end

% 统计结果
fprintf('总三角形数: %d\n', n_tris);
fprintf('正确外法线三角形数: %d (%.2f%%)\n',...
        sum(is_outward), 100*mean(is_outward));

% 可视化验证（显示前50个三角形）
figure;
hold on;
axis equal;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');

% 绘制球体节点
scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 10, 'b');

% 绘制三角形法线
% quiver3(centroids(1:50,1), centroids(1:50,2), centroids(1:50,3),...
%         normals(1:50,1), normals(1:50,2), normals(1:50,3),...
%         'AutoScale','on', 'Color','r');
quiver3(centroids(:,1), centroids(:,2), centroids(:,3),...
        normals(:,1), normals(:,2), normals(:,3),...
        'AutoScale','on', 'Color','r');

% 绘制错误法线（如果有）
if any(~is_outward)
    error_ids = find(~is_outward);
    quiver3(centroids(error_ids(1:min(5,end)),1),...
            centroids(error_ids(1:min(5,end)),2),...
            centroids(error_ids(1:min(5,end)),3),...
            normals(error_ids(1:min(5,end)),1),...
            normals(error_ids(1:min(5,end)),2),...
            normals(error_ids(1:min(5,end)),3),...
            'AutoScale','off', 'Color','k', 'LineWidth',2);
    legend({'节点','正确法线','错误法线'});
else
    legend({'节点','所有法线正确'});
end