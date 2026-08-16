%% ParticleDataAnalyze.m
%% 读取颗粒的受力和运动数据
clear
clc
%% 参数初始化
parameter;
DataPath = [Project_filefolder,'Data\',datafilefolder,'\OutputFile\'];

count = StepOutput;
time_steps = 0:count:EndStep; % 时间步序列
num_steps = length(time_steps);
step = StepSize;
nbody = num_body; % 根据实际情况设置刚体数量

%% 时间
time = zeros(num_steps,1);
%% 颗粒系统整体属性
Ek_move = zeros(num_steps,1); % 平动动能
Ek_spin = zeros(num_steps,1); % 转动动能
Ek = zeros(num_steps,1); %动能
Ep = zeros(num_steps,1); %势能

max_omg = zeros(num_steps,1);

%料斗卸料算例中颗粒残余
num_re = zeros(num_steps,1);
%% 某时刻所有颗粒的数据
%运动状态
pos = zeros(num_pt,3);
vel = zeros(num_pt,3);
quat = zeros(num_pt,4);
omg = zeros(num_pt,3);
omgsq = zeros(num_pt,3);
%受力数据
total_force = zeros(num_pt,3);
total_moment = zeros(num_pt,3);

FI_pt_data = zeros(num_pt, nbody*6);
%% 单个颗粒属性
k = 1000;
state = 1;
%运动状态
pos_ = zeros(num_steps,3);
vel_ = zeros(num_steps,3);
quat_ = zeros(num_steps,4);
omg_ = zeros(num_steps,3);
omgsq_ = zeros(num_steps,1);
%受力数据
total_force_ = zeros(num_steps,3);
total_moment_ = zeros(num_steps,3);
contact_force_ = zeros(num_steps,3);
contact_moment_ = zeros(num_steps,3);
wall_force_ = zeros(num_steps,3);
wall_moment_ = zeros(num_steps,3);

sd = 0;
%% 循环处理每个时间步
for i = 1:num_steps
    t = time_steps(i);
    time(i) = t*step;

    %% 读取粒子运动数据（ph文件）
    ph_file = fullfile(DataPath, ['state_particles\ph.', sprintf('%08d', t), '.bt']);
    tmp = load(ph_file);

    nb = size(tmp, 1);
    mass = tmp(:,3);
    radii = tmp(:,4);

    pos = tmp(:,5:7);
    vel = tmp(:,8:10);
    quat = tmp(:,11:14);
    omg = tmp(:,15:17);

    pos_(i,:) = tmp(k,5:7);
    vel_(i,:) = tmp(k,8:10);
    quat_(i,:) = tmp(k,11:14);
    omg_(i,:) = tmp(k,15:17);

    J=0.4*mass.*(radii.*radii);

    for j=1:nb
        max = omg(j,:)*omg(j,:)';
        if max>max_omg(i)
            max_omg(i) = max;
            sd = j;
        end
    end
    fprintf("max_id = %d\n",sd);
    sss = 0.0095;
    sss = 0.0195;
    sss = 0.0275;
    for j=1:nb
        Ek_move(i) = Ek_move(i) + 0.5*mass(j)*(vel(j,:)*vel(j,:)');
        Ek_spin(i) = Ek_spin(i) + 0.5*J(j)*(omg(j,:)*omg(j,:)');
        if pos(j,3)>0-sss
            num_re(i) = num_re(i) + 1;
        end
    end
    Ek(i)=Ek_move(i) + Ek_spin(i);



    %% 读取受力数据（F文件）
    % ph_file = fullfile(DataPath, ['force_particles\F.', sprintf('%08d', t), '.bt']);
    % F = load(ph_file);
    % 
    % total_force_(i,:) = F(k,16:18);
    % total_moment_(i,:) = F(k,19:21);
    % contact_force_(i,:) = F(k,1:3);
    % contact_moment_(i,:) = F(k,4:6);
    % wall_force_(i,:) = F(k,7:9);
    % wall_moment_(i,:) = F(k,10:12);


    
    % 合成总力和力矩
    % contact_force = [F{1}, F{2}, F{3}];
    % wall_force = [F{7}, F{8}, F{9}];
    % gravity = [F{13}, F{14}, F{15}];
    % total_force(i,:) = contact_force + wall_force + gravity;
    % 
    % contact_moment = [F{4}, F{5}, F{6}];
    % wall_moment = [F{10}, F{11}, F{12}];
    % total_moment(i,:) = contact_moment + wall_moment;

    % total_force(i,:) = [F{16}, F{17}, F{18}];
    % total_moment(i,:) = [F{19}, F{20}, F{21}];
    % % 提取新增的FI_pt数据
    % FI_pt_current = cell2mat(F(22:21+nbody*6))';
    % FI_pt_data(i,:) = FI_pt_current;
end

if state
    %% 绘图设置
    figure('Units','normalized','Position',[0.1 0.1 0.8 0.8])

    % 图1：速度分析
    subplot(2,2,1)
    plot(time, vecnorm(vel_,2,2), 'k-', 'LineWidth',2)
    hold on
    plot(time, vel_(:,1), 'r--')
    plot(time, vel_(:,2), 'g--')
    plot(time, vel_(:,3), 'b--')
    title(['粒子' num2str(k) '速度分析'])
    xlabel('时间(s)'), ylabel('速度 (m/s)')
    legend('合速度','X分量','Y分量','Z分量')
    grid on

    % 图2：位置轨迹
    subplot(2,2,2)
    hold on
    plot(time, pos_(:,1), 'r--')
    hold on
    plot(time, pos_(:,2), 'g--')
    plot(time, pos_(:,3), 'b--')
    title(['粒子' num2str(k) '位置分析'])
    xlabel('时间(s)'), ylabel('位置 (m)')
    legend('X分量','Y分量','Z分量')
    grid on

    % 图3：角速度分析
    subplot(2,2,3)
    plot(time, vecnorm(omg_,2,2), 'k-', 'LineWidth',2)
    hold on
    plot(time, omg_(:,1), 'r--')
    plot(time, omg_(:,2), 'g--')
    plot(time, omg_(:,3), 'b--')
    title(['粒子' num2str(k) '角速度分析'])
    xlabel('时间(s)'), ylabel('速度 (rad/s)')
    legend('合角速度','X分量','Y分量','Z分量')
    grid on

    % 图4：三维运动轨迹
    subplot(2,2,4)
    plot3(pos_(:,1), pos_(:,2), pos_(:,3), 'm-')
    title(['粒子' num2str(k) '运动轨迹'])
    xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
    axis([0, 0.5, 0, 0.5, 0 , 0.5]);
    axis manual; 
    %axis equal
    grid on

    figure('Units','normalized','Position',[0.1 0.1 0.8 0.8])
    % 图1：受力分析
    subplot(2,2,1)
    force_magnitude = vecnorm(total_force_,2,2);
    plot(time, force_magnitude, 'k-', 'LineWidth',2)
    hold on
    plot(time, total_force_(:,1), 'r--')
    plot(time, total_force_(:,2), 'g--')
    plot(time, total_force_(:,3), 'b--')
    title(['粒子' num2str(k) '受力分析'])
    xlabel('时间(s)'), ylabel('力 (N)')
    legend('合力','X分量','Y分量','Z分量')
    grid on

    % 图2：力矩分析
    subplot(2,2,2)
    moment_magnitude = vecnorm(total_moment_,2,2);
    plot(time, moment_magnitude, 'k-', 'LineWidth',2)
    hold on
    plot(time, total_moment_(:,1), 'r--')
    plot(time, total_moment_(:,2), 'g--')
    plot(time, total_moment_(:,3), 'b--')
    title(['粒子' num2str(k) '力矩分析'])
    xlabel('时间(s)'), ylabel('力矩 (N·m)')
    legend('合力矩','X分量','Y分量','Z分量')
    grid on

    % 图1：受力分析
    subplot(2,2,3)
    force_magnitude = vecnorm(wall_force_,2,2);
    plot(time, force_magnitude, 'k-', 'LineWidth',2)
    hold on
    plot(time, wall_force_(:,1), 'r--')
    plot(time, wall_force_(:,2), 'g--')
    plot(time, wall_force_(:,3), 'b--')
    title(['粒子' num2str(k) '受力分析-墙壁'])
    xlabel('时间(s)'), ylabel('力 (N)')
    legend('合力','X分量','Y分量','Z分量')
    grid on

    % 图4：力矩分析
    subplot(2,2,4)
    moment_magnitude = vecnorm(wall_moment_,2,2);
    plot(time, moment_magnitude, 'k-', 'LineWidth',2)
    hold on
    plot(time, wall_moment_(:,1), 'r--')
    plot(time, wall_moment_(:,2), 'g--')
    plot(time, wall_moment_(:,3), 'b--')
    title(['粒子' num2str(k) '墙壁接触力矩分析'])
    xlabel('时间(s)'), ylabel('力矩 (N·m)')
    legend('合力矩','X分量','Y分量','Z分量')
    grid on

    figure('Units','normalized','Position',[0.1 0.1 0.8 0.8])
    % 图1：受力分析
    subplot(2,2,1)
    force_magnitude = vecnorm(contact_force_,2,2);
    plot(time, force_magnitude, 'k-', 'LineWidth',2)
    hold on
    plot(time, contact_force_(:,1), 'r--')
    plot(time, contact_force_(:,2), 'g--')
    plot(time, contact_force_(:,3), 'b--')
    title(['粒子' num2str(k) '受力分析-接触p'])
    xlabel('时间(s)'), ylabel('力 (N)')
    legend('合力','X分量','Y分量','Z分量')
    grid on

    % 图2：力矩分析
    subplot(2,2,2)
    moment_magnitude = vecnorm(contact_moment_,2,2);
    plot(time, moment_magnitude, 'k-', 'LineWidth',2)
    hold on
    plot(time, contact_moment_(:,1), 'r--')
    plot(time, contact_moment_(:,2), 'g--')
    plot(time, contact_moment_(:,3), 'b--')
    title(['粒子' num2str(k) '力矩分析-接触p'])
    xlabel('时间(s)'), ylabel('力矩 (N·m)')
    legend('合力矩','X分量','Y分量','Z分量')
    grid on


end

figure
plot(time,max_omg)

figure
hold on
plot(time, Ek, 'k-', 'LineWidth',2)
plot(time, Ek_move, 'b-', 'LineWidth',2)
plot(time, Ek_spin, 'r-', 'LineWidth',2)


a = load("experiment.csv");
b = load("simulation.csv");
c = num_re/num_re(1);
d = abs(a-c);

figure('Position', [100, 100, 490, 350]); 
hold on; box on;%grid on;
plot(time, c, 'k-')
plot(a(:,1), a(:,2), 'b-')
plot(b(:,1), b(:,2), 'r-')
plot()
xticks([0, 2, 4, 6, 8]);% 设置X轴刻度仅在 0, 10, 20, 30, 40, 50 处显示
set(gca, 'XGrid', 'off', 'YGrid', 'off');
xlabel('Time (s)');   % \mus 表示微秒符号
ylabel('Residual particle fraction'); % \mum 表示微米符号
lgd = legend('Sim-ours','Experiment','Sim-cited','Location', 'northeast');
%set(lgd, 'FontSize', 14);
% 统一设置字体为 Times New Roman，字号为 12
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);