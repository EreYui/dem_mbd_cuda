% parameter

%% 积分器
StepOutput = 500;       % 输出文件步长间隔
StartStep  = 0;         % 初始步
EndStep	   = 100000;    % 结束步
StepSize   = 2.5e-5;   % 步长

%% 颗粒
num_pt = 25215;                        % 颗粒数量
file_num_pt = EndStep/StepOutput + 1;  % 颗粒文件数量
%% 刚体
num_body = 1;                          % 刚体数量
file_num_body = EndStep/StepOutput;    % 刚体文件数量

%% 文件
Project_filefolder = 'F:\paper1code\dem_mbd_omp\';        % 项目文件夹
datafilefolder     = 'DATA5';               % 本算例输入输出文件夹
picture_num        = EndStep/StepOutput;   % 输出图片数量

%% 物理参数
g = 9.8;

