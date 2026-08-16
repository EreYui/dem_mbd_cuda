%% 将paraview输出的图片组合为视频
clear
clc

Project_filefolder = 'F:\paper1code\dem_mbd_omp\';        % 项目文件夹
datafilefolder     = 'DATA5';               % 本算例输入输出文件夹
picture_num        = 641;   % 输出图片数量

framesPath = [Project_filefolder,'Data\',datafilefolder,'\OutputFile\pic\'];%图像序列所在路径，同时要保证图像大小相同
videoName  = [Project_filefolder,'Data\',datafilefolder,'\OutputFile\video001.avi'];%表示将要创建的视频文件的名字
fps = 20; %帧率
startFrame = 0; %从哪一帧开始
endFrame   = picture_num; %哪一帧结束
 
if(exist('videoName','file'))
    delete videoName.avi
end
 
%生成视频的参数设定
aviobj=VideoWriter(videoName);  %创建一个avi视频文件对象，开始时其为空
aviobj.FrameRate=fps;
open(aviobj);%Open file for writing video data
%读入图片
for i=startFrame:2:endFrame
    fileName=sprintf('pic.%04d',i);    %根据文件名而定 我这里文件名是pic.0000.png pic.0001.png pic.0002.png ....
    frames=imread([framesPath,fileName,'.png']);
    writeVideo(aviobj,frames);
    fprintf('正在处理: %d/%d\n', i, picture_num);
end
close(aviobj);% 关闭创建视频
disp('视频已输出。');