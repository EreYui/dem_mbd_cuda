%% 将输出的ph.xxx.bt格式转换为vtk格式
%% 这里只输出一个文件，用于查看生成的颗粒文件
clc
clear
%%
version='# vtk DataFile Version 3.0\n';
output='ps-dem information visualization\n';
Ascii='ASCII\n';
Dataset='DATASET UNSTRUCTURED_GRID\n';
Velocity='VECTORS velocity double';
Spin='VECTORS spin double';

FILENAME = 'generate1';
%
resfile = [FILENAME,'.txt'];
%
tmp = load(resfile);
mass = tmp(:,3);
radii = tmp(:,4);
pos = tmp(:,5:7);
vel = tmp(:,8:10);
quat = tmp(:,11:14);
omg = tmp(:,15:17);
type=tmp(:,2);
nb = length(tmp(:,1));
input= [FILENAME,'.vtk'];
pointNum=['POINTS ',num2str(nb),' double\n'];
cellNum=['CELLS ',num2str(nb),' ',num2str(2*nb),'\n'];
cellType=['CELL_TYPES ',num2str(nb),'\n'];
point_Data=['POINT_DATA  ',num2str(nb),'\n'];
fieldData=['FIELD FieldData 4 ',num2str(nb),'\n'];
r=['radius 1 ',num2str(nb),' double\n'];
color=['color 1 ',num2str(nb),' int\n'];
%
J=0.4*mass.*(radii.*radii);

fid = fopen(input,'wt+');
fprintf(fid,version);
fprintf(fid,output);
fprintf(fid,Ascii);
fprintf(fid,Dataset);
fprintf(fid,pointNum);
Formats='%22.15e %22.15e %22.15e \n';
for j=1:nb
    fprintf(fid,Formats,pos(j,1),pos(j,2),pos(j,3));
end
fprintf(fid,cellNum);
Formats='%5d %5d \n';
for j=1:nb
    fprintf(fid,Formats,1,j-1);
end
fprintf(fid,cellType);
Formats='%5d \n';
for j=1:nb
    fprintf(fid,Formats,1);
end
%%
fprintf(fid,'POINT_DATA %d\n',nb);
fprintf(fid,'FIELD FieldData 4\n');
fprintf(fid,'number 1 %d int\n',nb);
for ii = 1:nb
    fprintf(fid,'%d\n',ii-1);
end
%%
fprintf(fid,r);
for ii = 1:nb
    fprintf(fid,'%f\n',radii(ii));
end
%%
fprintf(fid,'velocity 1 %d double\n',nb);
for ii = 1:nb
    fprintf(fid,'%f\n',norm(vel(ii,:),2));
end
fprintf(fid,'angular_velocity 1 %d double\n',nb);
for ii = 1:nb
    fprintf(fid,'%f\n',norm(omg(ii,:),2));
end
%%

fclose(fid);
%%
fclose all;
