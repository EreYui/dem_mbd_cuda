%% 将输出的ph.xxx.bt格式转换为vtk格式
%% 处理颗粒的输出数据
clc
clear

parameter;
%%
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

% 预定义VTK文件头信息
version = '# vtk DataFile Version 3.0\n';
output = 'ps-dem information visualization\n';
Ascii = 'ASCII\n';
Dataset = 'DATASET UNSTRUCTURED_GRID\n';

%%
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

    % 生成输出文件名
    vtkfile = fullfile(vtkfolder, ['ph.', sprintf('%08d', ifile), '.vtk']);
    fid = fopen(vtkfile, 'wt+');



    % 写入VTK头信息
    fprintf(fid, version);
    fprintf(fid, output);
    fprintf(fid, Ascii);
    fprintf(fid, Dataset);
    
    % 写入点坐标（向量化操作）
    fprintf(fid, 'POINTS %d double\n', nb);
    fprintf(fid, '%22.15e %22.15e %22.15e\n', pos'); % 转置矩阵加速写入
    
    % 写入单元信息（向量化操作）
    fprintf(fid, 'CELLS %d %d\n', nb, 2*nb);
    fprintf(fid, '%d %d\n', [ones(nb,1), (0:nb-1)']'); % 合并数据加速写入
    
    fprintf(fid, 'CELL_TYPES %d\n', nb);
    fprintf(fid, '%d\n', ones(nb,1)); % 所有单元类型为1
    
    % 写入点数据（向量化操作）
    fprintf(fid, 'POINT_DATA %d\n', nb);
    
    % 编号字段
    fprintf(fid, 'FIELD FieldData 4\n');
    fprintf(fid, 'number 1 %d int\n', nb);
    fprintf(fid, '%d\n', (0:nb-1)');
    
    % 半径字段
    fprintf(fid, 'radius 1 %d double\n', nb);
    fprintf(fid, '%f\n', radii);
    
    % 速度模长字段
    fprintf(fid, 'velocity 1 %d double\n', nb);
    fprintf(fid, '%f\n', sqrt(sum(vel.^2, 2)));
    
    % 角速度模长字段
    fprintf(fid, 'angular_velocity 1 %d double\n', nb);
    fprintf(fid, '%f\n', sqrt(sum(omg.^2, 2)));
    
    fclose(fid);
    fprintf('Processed file: %d/%d\n', idx, numel(fileIndices));
end

%%
fclose all;

fig = figure('Name',[filefolder,'\Kinetic Energy']);
box on;
hold on
plot(N*step,Ek,'k','LineWidth',2);
plot(N*step,Ek_spin,'r','LineWidth',2);
plot(N*step,Ek_move,'b','LineWidth',2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Time(s)','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
ylabel('Kinetic Energy','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
legend('Ek', 'Ek_spin', 'Ek_move', 'Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex')
grid on;
print(gcf, [filefolder,'\Kinetic Energy'], '-dpng', '-r300'); % 保存为300 DPI的PNG格式
% 保存成矢量图（SVG、PDF）
saveas(fig, fig.Name, 'svg');
saveas(fig, fig.Name, 'pdf');
