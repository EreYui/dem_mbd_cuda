
% 合并两个数据文件，重新编号并保持其他数据不变
% 输入:
%   file1 - 第一个输入文件名 (字符串)
%   file2 - 第二个输入文件名 (字符串)
%   output_file - 输出文件名 (字符串)

file1 = "pt1.txt";
file2 = "pt2.txt";
output_file = "output.txt";

    % 读取数据
    A = readmatrix(file1);
    B = readmatrix(file2);

    % 校验列数
    if size(A,2) ~= 17 || size(B,2) ~= 17
        error('数据必须为17列');
    end

    % 重新编号
    nA = size(A,1);
    B(:,1) = B(:,1) + nA;   % 接在A之后

    % 合并
    C = [A; B];

    % 写入文件（空格分隔，无逗号）
    fid = fopen(output_file, 'w');
    for i = 1:size(C,1)
        fprintf(fid, '%g', C(i,1));
        for j = 2:17
            fprintf(fid, ' %g', C(i,j));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);

fprintf('合并完成！\n');
% fprintf('文件1: %s (%d 行)\n', file1, size(data1, 1));
% fprintf('文件2: %s (%d 行)\n', file2, size(data2, 1));
% fprintf('合并后: %s (%d 行)\n', output_file, size(merged_data, 1));
% fprintf('编号范围: %d 到 %d\n', merged_data(1,1), merged_data(end,1));
