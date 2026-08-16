#include "wall.h"
#include "matrixtrans.h"
#include "file_utils.h"
#include <cassert>
#include <sstream>
#include <fstream>

void WALL::LoadWalls(const std::string& parafile) {
    const auto wallPath = file_utils::requireInputFile(
        parafile, "wall input file / 墙体输入文件");
    FILE* infile = file_utils::openInputFile(
        wallPath.u8string(), "wall input file / 墙体输入文件");
    if (fscanf(infile, "%d", &number) != 1 || number < 0) {
        fclose(infile);
        throw std::runtime_error(
            "Invalid wall count / 墙体数量格式错误: "
            + file_utils::absolutePathForMessage(wallPath));
    }
    const int Number = number;

    // Alloc wall
    Orig = new vector3d[Number];
    N = new vector3d[Number];
    for (int j = 0; j < Number; j++) {
        if (fscanf(infile, "%lf%lf%lf%lf%lf%lf",
                &Orig[j][0], &Orig[j][1], &Orig[j][2],
                &N[j][0], &N[j][1], &N[j][2]) != 6) {
            fclose(infile);
            throw std::runtime_error(
                "Invalid wall record / 墙体数据格式错误: "
                + file_utils::absolutePathForMessage(wallPath)
                + ", record / 记录 "
                + std::to_string(j + 1));
        }
    }

    fclose(infile);
}
