%% 处理计算完成后刚体系的数据，用于分析
clear
clc
parameter;
%% 用户输入参数
count = StepOutput;
n = file_num_body;       % 文件数量
nb = num_body;        % 刚体数量
nd = n;  % 数据行数

mass = 3;

ptfolder = [Project_filefolder,'Data\',datafilefolder,'\OutputFile\state_bodys\'];
% 预分配数据存储数组
data = zeros(nb, nd , 20);
%% 主逻辑
for timeStep = 1:n
    % 生成文件名
    resfile = fullfile(ptfolder, ['At.', sprintf('%08d', timeStep * count), '.bt']);
    tmp = load(resfile);
    
    % 处理每个刚体数据
    for bodyID = 1:nb
        row_num = bodyID;  % 计算目标行号
        line_data = tmp(bodyID,:); % 解析行数据
        
        % 提取所需数据
        data(bodyID,timeStep, :) = line_data;
    end
end

k = 1;
%% 位置
figure
box on;
hold on
plot(data(k,:,1),data(k,:,2),'k','LineWidth',2);
plot(data(k,:,1),data(k,:,3),'r','LineWidth',2);
plot(data(k,:,1),data(k,:,4),'b','LineWidth',2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Time(s)','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
ylabel('Position','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
grid on;
legend('x', 'y', 'z', 'Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');
%% 速度
figure
box on;
hold on
plot(data(k,:,1),data(k,:,5),'k','LineWidth',2);
plot(data(k,:,1),data(k,:,6),'r','LineWidth',2);
plot(data(k,:,1),data(k,:,7),'b','LineWidth',2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Time(s)','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
ylabel('Position','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
grid on;
legend('vx', 'vy', 'vz', 'Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');
%% 角速度
figure
box on;
hold on
plot(data(k,:,1),data(k,:,8),'k','LineWidth',2);
plot(data(k,:,1),data(k,:,9),'r','LineWidth',2);
plot(data(k,:,1),data(k,:,10),'b','LineWidth',2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Time(s)','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
ylabel('AngularVelocity','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
grid on;
legend('omgx', 'omgy', 'omgz', 'Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');
%% shouli
figure
box on;
hold on
plot(data(k,:,1),data(k,:,15),'k','LineWidth',2);
plot(data(k,:,1),data(k,:,16),'r','LineWidth',2);
plot(data(k,:,1),data(k,:,17),'b','LineWidth',2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Time(s)','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
ylabel('Force','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
grid on;
legend('Fx', 'Fy', 'Fz', 'Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');

%% 力矩
figure
box on;
hold on
plot(data(k,:,1),data(k,:,18),'k--','LineWidth',2);
plot(data(k,:,1),data(k,:,19),'r--','LineWidth',2);
plot(data(k,:,1),data(k,:,20),'b--','LineWidth',2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Time(s)','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
ylabel('Force','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
grid on;
legend('Tx', 'Ty', 'Tz', 'Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');

%% 能量
Ek = 0.5*mass*(data(k,:,5).*data(k,:,5) + data(k,:,6).*data(k,:,6) + data(k,:,7).*data(k,:,7))...
    +0.5*3e-4*(data(k,:,8).*data(k,:,8) + data(k,:,9).*data(k,:,9) + data(k,:,10).*data(k,:,10));
Ep = mass * 9.8 * data(k,:,4);
E = Ek + Ep;
figure
box on;
hold on
plot(data(k,:,1),Ek,'k','LineWidth',2);
plot(data(k,:,1),Ep,'r','LineWidth',2);
plot(data(k,:,1),E,'b','LineWidth',2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Time(s)','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
ylabel('Force','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
grid on;
legend('Ek', 'Ep', 'E', 'Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');


%% 
% en = zeros(n,1);
% for i = 1:n
%     en(i) = 0.5*m*(data(1,i,5:7).norm)^
% end
% %% 创建结果表格并设置类型
% resultTable = array2table(data,...
%     'VariableNames', {'time','x','y','z','vx','vy','vz','ax','ay','az','q1','q2','q3','q4','f1','f2','f3','f4','f5','f6'});
% % resultTable.timeStep = int32(resultTable.timeStep);
% % resultTable.bodyID = int32(resultTable.bodyID);
% 
% %% 写入CSV文件
% writetable(resultTable, 'output_data.csv');
figure%shouli
box on;
hold on
plot(data(k,1:64,1),data(k,1:64,15),'k','LineWidth',2);
plot(data(k,1:64,1),data(k,1:64,16),'r','LineWidth',2);
plot(data(k,1:64,1),data(k,1:64,17),'b','LineWidth',2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Time(s)','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
ylabel('Force','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
grid on;
legend('Fx', 'Fy', 'Fz', 'Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');

%% 位置
figure
box on;
hold on
plot(data(k,1:64,1),data(k,1:64,2),'k','LineWidth',2);
plot(data(k,1:64,1),data(k,1:64,3),'r','LineWidth',2);
plot(data(k,1:64,1),data(k,1:64,4),'b','LineWidth',2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('Time(s)','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
ylabel('Position','FontSize', 14, 'FontName', 'Times New Roman','FontWeight','bold')
grid on;
legend('x', 'y', 'z', 'Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');