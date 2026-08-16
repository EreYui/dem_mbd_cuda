clear
clc
%% 
% 按照给定空间生成颗粒
% 颗粒之间互相不接触
% 需要做一下沉降
output_filename = 'pt1.txt'; % 输出文件名
%% 
x = [0 1];
y = [0 1];
z = [0 0.6];
d = 0.02;% 颗粒直径

% x = [-0.05 0.05];
% y = [-0.05 0.05];
% z = [-0.0195 0.0275];
% d = 0.02;% 颗粒直径

%z = z+0.1;

% 每个方向的长度
wx = x(2)-x(1);
wy = y(2)-y(1);
wz = z(2)-z(1);

ka = 1.07;

% 将空间划分为边长为d*ka的网格，得到的网格数量
xnum = floor(wx/(d*ka));
ynum = floor(wy/(d*ka));
znum = floor(wz/(d*ka));

kax = wx / xnum / d;
kay = wy / ynum / d;
kaz = wz / znum / d;

den = 1357;

nde = 59248-57694;%需要随机删除的颗粒数量
n1 = 0;%已删除的颗粒数量
% tag1 = 0;

n = xnum*znum*ynum;
n = 57694;
data = zeros(n, 17);
ib = 0;

for i = 1:xnum
    for j = 1:ynum
        for k = 1:znum
            if n1<nde
                ttt = rand;
                if ttt<0.05
                    n1 = n1 + 1;
                    continue;
                end
            else
                tag = 1;
            end
            ib = ib + 1;
            %ib = (k-1)*xnum*ynum + (j-1)*xnum + i;
            data(ib,1) = ib-1;%num
            data(ib,2) =0;%state
            data(ib, 3) = pi*d^3*den/6;%mass
            data(ib, 4) = d/2;%r
            data(ib, 5) = x(1) + (i-1)*d*kax + d/2 + (kax - 1)*d*rand(1)/2;%position,x,y,z
            data(ib, 6) = y(1) + (j-1)*d*kay + d/2 + (kay - 1)*d*rand(1)/2;
            data(ib, 7) = z(1) + (k-1)*d*kaz + d/2 + (kaz - 1)*d*rand(1)/2;
            data(ib, 8:10) = zeros(1,3);%velocity,vx,vy,vz
            data(ib, 11) = 1;
            data(ib, 12:14) = zeros(1,3);%四元数
            data(ib, 15:17) = zeros(1,3);%角速度
        end
    end
end

%% 边界接触检测
% 参数设置
x_wall = x;  % x方向墙坐标
y_wall = y;  % y方向墙坐标
r = d / 2;         % 颗粒半径


x = data(:, 5);  % 生成[-1.5,1.5]随机x坐标
y = data(:, 6);    % 生成[-2,2]随机y坐标
num_particles = size(x,1);
% 向量化接触判断
contact = x < (x_wall(1) + r) | ...   % 左墙接触
          x > (x_wall(2) - r) | ...   % 右墙接触
          y < (y_wall(1) + r) | ...   % 下墙接触
          y > (y_wall(2) - r);        % 上墙接触

% 结果分析
if any(contact)
    fprintf('发现边界接触! 受影响颗粒数量: %d (%.2f%%)\n',...
            sum(contact), sum(contact)/num_particles*100);
    
    % 显示首个接触颗粒信息
    first_idx = find(contact,1);
    fprintf('首接触颗粒信息: x=%.4f, y=%.4f (半径=%.4f)\n',...
            x(first_idx), y(first_idx), r);
    
    % 可视化接触分布（可选）
    figure
    scatter(x(contact), y(contact), 10, 'r','filled')
    hold on
    rectangle('Position',[x_wall(1) y_wall(1) diff(x_wall) diff(y_wall)],...
              'EdgeColor','b','LineWidth',2)
    title('颗粒接触边界分布图')
    xlabel('X'), ylabel('Y')
    axis equal
else
    disp('所有颗粒均未接触边界')
end

%% 颗粒接触检测

% 合并坐标矩阵
points = data(:, 5:7);
total_particles = size(points, 1);

% 构建KD树加速邻居搜索
kdtree = KDTreeSearcher(points);

% 分块处理参数设置
block_size = 10000;  % 根据内存调整块大小
has_contact = false;

% 分块遍历所有颗粒
for start_idx = 1:block_size:total_particles
    end_idx = min(start_idx+block_size-1, total_particles);
    current_indices = start_idx:end_idx;
    
    % 批量查询当前块内所有点的邻居
    neighbor_cells = rangesearch(kdtree, points(current_indices,:), d);
    
    % 检查每个点的邻居列表
    for i = 1:length(neighbor_cells)
        global_idx = current_indices(i);
        neighbors = neighbor_cells{i};
        
        % 移除自身并检查实际距离
        neighbors(neighbors == global_idx) = [];
        if ~isempty(neighbors)
            % 精确计算距离（避免网格边缘误判）
            distances = sqrt(sum((points(neighbors,:) - points(global_idx,:)).^2, 2));
            if any(distances < d)
                has_contact = true;
                break;
            end
        end
    end
    
    if has_contact
        fprintf('发现接触! 首对接触颗粒：%d 与 %d\n', global_idx, neighbors(find(distances<d,1)));
        break;
    end
end

% 输出最终结果
if has_contact
    disp('=== 存在接触的颗粒 ===');
else
    disp('=== 所有颗粒均未接触 ===');
end

%% 保存新数据到文件
writematrix(data, output_filename, 'Delimiter', ' '); % 以空格分隔保存

disp(['数据处理完成并已保存到 ',output_filename]);


   