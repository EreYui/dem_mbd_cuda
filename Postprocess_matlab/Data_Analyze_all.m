% Data_Analyze_all
%% 处理计算完成后颗粒和刚体系的数据，用于分析
clear
clc
parameter;
%% 用户输入参数


%% 获取刚体（系）数据
nd = file_num_body;  % 数据行数
mass = 0.1;

ptfolder = [Project_filefolder,'Data\',datafilefolder,'\OutputFile\state_bodys\'];
% 预分配数据存储数组
data_bodys = zeros(num_body, nd , 20);
%% 主逻辑,提取刚体的数据
for timeStep = 1:file_num_body
    % 生成文件名
    resfile = fullfile(ptfolder, ['At.', sprintf('%08d', timeStep * StepOutput), '.bt']);
    tmp = load(resfile);
    
    % 处理每个刚体数据
    for bodyID = 1:num_body
        row_num = bodyID;  % 计算目标行号
        line_data = tmp(bodyID,:); % 解析行数据
        
        % 提取所需数据
        data_bodys(bodyID,timeStep, :) = line_data;
    end
end
k = 1;
Ek_body = 0.5*mass*(data_bodys(k,:,5).*data_bodys(k,:,5) + data_bodys(k,:,6).*data_bodys(k,:,6) + data_bodys(k,:,7).*data_bodys(k,:,7))...
    +0.5*3e-4*(data_bodys(k,:,8).*data_bodys(k,:,8) + data_bodys(k,:,9).*data_bodys(k,:,9) + data_bodys(k,:,10).*data_bodys(k,:,10));
Ep_body = mass * 9.8 * data_bodys(k,:,4);
E_body = Ek_body + Ep_body;
E_body = [E_body(1),E_body];
%% 提取颗粒的数据
step = StepSize;
fileseries = [0, StepOutput, EndStep];
fileIndices = fileseries(1):fileseries(2):fileseries(3); % 生成所有需要处理的文件索引
%%
filefolder = [Project_filefolder,'Data\',datafilefolder,'\OutputFile\'];
ptfolder = [filefolder,'state_particles\'];
vtkfolder = [filefolder,'vtk_particles\'];

N=fileseries(1):fileseries(2):fileseries(3);
Ek=zeros(length(N),1); 
Ek_move=zeros(length(N),1); 
Ek_spin=zeros(length(N),1); 
E = zeros(length(N),1);
% 创建并行池（根据系统资源调整池大小）
if isempty(gcp('nocreate'))
    parpool;
end

parfor idx = 1:numel(fileIndices)
    ifile = fileIndices(idx);
    resfile = fullfile(ptfolder, ['ph.', sprintf('%08d', ifile), '.bt']);
    tmp = load(resfile);
    
    % 提取数据
    nb = size(tmp, 1);
    mass = tmp(:,3);
    radii = tmp(:,4);
    pos = tmp(:,5:7);
    vel = tmp(:,8:10);
    omg = tmp(:,15:17);

    J=0.4*mass.*(radii.*radii);

    for j=1:nb
        Ek_move(idx) = Ek_move(idx) + 0.5*mass(j)*(vel(j,:)*vel(j,:)');
        Ek_spin(idx) = Ek_spin(idx) + 0.5*J(j)*(omg(j,:)*omg(j,:)');
        
    end
    Ek(idx)=Ek_move(idx) + Ek_spin(idx);
    E(idx) = Ek(idx) + sum(mass.*pos(:,3)*g);
end

%% 绘制系统能量图
fig = figure('Name',[filefolder,'\Kinetic Energy']);
box on;
hold on
plot(N*step,E,'k','LineWidth',2);
plot(N*step,Ek_spin,'r','LineWidth',2);
plot(N*step,Ek_move,'b','LineWidth',2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Time(s)','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
ylabel('Kinetic Energy','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
legend('Ek', 'Ek_spin', 'Ek_move', 'Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex')
grid on;

all_E = E + E_body;
fig = figure('Name',[filefolder,'\all Kinetic Energy']);
box on;
hold on
plot(N*step,E,'k','LineWidth',2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Time(s)','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
ylabel('Kinetic Energy','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
legend('Ek', 'Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex')
grid on;