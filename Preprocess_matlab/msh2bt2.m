%% 将Gmsh导出的msh文件转化为bt文件，方便仿真程序读取

clear
name = 'cube_150mm_surface';
convertMSHtoBT([name,'.msh'],[name,'.bt']);

function convertMSHtoBT(mshFileName, btFileName)
    % 打开MSH文件
    fileID = fopen(mshFileName, 'r');
    if fileID == -1
        error('无法打开MSH文件: %s', mshFileName);
    end
    
    % 初始化变量
    nodes = [];
    triangles = [];
    readingNodes = false;
    readingElements = false;
    readingTriangles = false;
    
    % 逐行读取MSH文件
    tline = fgetl(fileID);
    while ischar(tline)
        % 去除行尾的换行符和空格
        tline = strtrim(tline);
        
        % 检查是否开始读取节点
        if strcmp(tline, '$Nodes')
            readingNodes = true;
            readingTriangles = false; % 确保不在同时读取三角形
            tline = fgetl(fileID);
            continue;
        end
        
        if strcmp(tline, '$Elements')
            readingNodes = false;
            readingElements = true; % 确保不在同时读取三角形
            tline = fgetl(fileID);
            continue;
        end
        
        % 如果在读取节点，则处理节点数据
        if readingNodes
            % 检查行是否只包含一个整数（节点ID）或三个双精度数（节点坐标）
            
            if (numel(strsplit(tline)) == 3) % 检查是否有空格，假设坐标行包含空格分隔的值
                data = strsplit(tline);
                % 检查数据是否包含四个元素（可能是ID和三个坐标）
                    % 提取坐标（假设最后三个是双精度数）
                    x = str2double(data{1});
                    y = str2double(data{2});
                    z = str2double(data{3});
                    % 节点ID可能不是必需的，或者可以在之前的行中读取并存储在一个数组中
                    % 如果需要ID，则可以从data{1}中提取
                    nodes = [nodes; x, y, z]; % 只存储坐标
             end
        end
        
        if readingElements
            % 检查行是否只包含4个整数（节点ID）或三个节点编号
            data = strsplit(tline);
             if (str2double(data{1}) == 2) && (numel(strsplit(tline)) == 4) % 单个整数（可能是节点ID，但通常不是坐标行）
                readingTriangles = true ;
                tline = fgetl(fileID);
                continue; 
             end
            if  readingTriangles && (numel(strsplit(tline)) == 4) % 单个整数（可能是节点ID，但通常不是坐标行）
                    n1 = str2double(data{2});
                    n2 = str2double(data{3});
                    n3 = str2double(data{4});
                    triangles = [triangles; n1, n2, n3]; % 只存储坐标
             end
        end
        
        
        if strcmp(tline, '$EndElements')
            break;
        end
        
        % 读取下一行
        tline = fgetl(fileID);

    end
    
    fclose(fileID);
    
[nP,~]=size(nodes(:,1));[nF,~]=size(triangles(:,1));
file=fopen(btFileName,'w');
if file == -1
    error('无法打开文件'); % 错误处理
end

Px=zeros(nP,1);Py=zeros(nP,1);Pz=zeros(nP,1);
for i=1:nP
    Px(i)=nodes(i,1);
    Py(i)=nodes(i,2);
    Pz(i)=nodes(i,3);
end
%%
fprintf(file,'%d\t%d\n',nP,nF);
for i=1:nP
    fprintf(file,'%.15e\t%.15e\t%.15e\n',Px(i)/1000,Py(i)/1000,Pz(i)/1000);
end
for i=1:nF
    fprintf(file,'%d\t%d\t%d\n',triangles(i,1),triangles(i,2),triangles(i,3));
end
fclose(file);
    
    
end