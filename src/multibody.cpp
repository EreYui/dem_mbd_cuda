#include "multibody.h"
#include "constant.h"
#include "file_utils.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cctype>

using namespace std;

inline bool LoadBodysFromCSV(BODYSET& bodyset, const std::string& filename);

void BODYSET::LoadBodys(const std::string& filename)
{
    BODYSET& bodyset = *this;
    if (LoadBodysFromCSV(*this, filename)) {
        std::cout << " ɹ " << Num << "" << std::endl;

        // 输出刚体信息
        for (int i = 0; i < Num; i++) {
            std::cout << "\n刚体 #" << i + 1 << " - " << body[i].name;
            std::cout << "\n 状态 state: " << body[i].state;
            std::cout << "\n 质量: " << body[i].Mass;
            std::cout << "\n 质心位置: ("
                << body[i].MassCenter[0] << ", "
                << body[i].MassCenter[1] << ", "
                << body[i].MassCenter[2] << ")";
            std::cout << "\n 姿态 (四元数): ("
                << body[i].orien[0] << ", "
                << body[i].orien[1] << ", "
                << body[i].orien[2] << ", "
                << body[i].orien[3] << ")";
            std::cout << "\n 惯性矩阵:\n" << body[i].I << std::endl;
        }
    }
    else {
        throw std::runtime_error(
            "Failed to load rigid bodies / 刚体列表加载失败: " + filename);
    }

    for (int i = 0; i < bodyset.Num; i++)
    {
        double radius = 0.0;
        vector3d I{ 0.0, 0.0, 0.0 };
        const auto meshPath = file_utils::requireInputFile(
            bodyset.body[i].name,
            "rigid-body mesh file / 刚体网格文件 #"
                + std::to_string(i + 1));
        FILE* infile = file_utils::openInputFile(
            meshPath.u8string(),
            "rigid-body mesh file / 刚体网格文件 #"
                + std::to_string(i + 1));

        // read vertex and triangle counts
        if (fscanf(infile, "%d%d", &bodyset.body[i].pointNum,
                &bodyset.body[i].Size) != 2
            || bodyset.body[i].pointNum <= 0
            || bodyset.body[i].Size <= 0) {
            fclose(infile);
            throw std::runtime_error(
                "Invalid mesh header / 刚体网格头格式错误: "
                + file_utils::absolutePathForMessage(meshPath));
        }
        
        // allocate vertex/triangle buffers
        bodyset.body[i].point = new vector3d[bodyset.body[i].pointNum];
        bodyset.body[i].tri = new vector3i[bodyset.body[i].Size];
        // read vertices
        for (int k = 0; k < bodyset.body[i].pointNum; ++k) {
            if (fscanf(infile, "%lf%lf%lf",
                    &bodyset.body[i].point[k][0],
                    &bodyset.body[i].point[k][1],
                    &bodyset.body[i].point[k][2]) != 3) {
                fclose(infile);
                throw std::runtime_error(
                    "Invalid mesh vertex / 刚体网格顶点格式错误: "
                    + file_utils::absolutePathForMessage(meshPath)
                    + ", vertex / 顶点 "
                    + std::to_string(k + 1));
            }
        }
        // read triangles
        for (int k = 0; k < bodyset.body[i].Size; ++k) {
            if (fscanf(infile, "%d%d%d",
                    &bodyset.body[i].tri[k][0],
                    &bodyset.body[i].tri[k][1],
                    &bodyset.body[i].tri[k][2]) != 3) {
                fclose(infile);
                throw std::runtime_error(
                    "Invalid mesh face / 刚体网格面格式错误: "
                    + file_utils::absolutePathForMessage(meshPath)
                    + ", face / 面 "
                    + std::to_string(k + 1));
            }
            for (int vertex = 0; vertex < 3; ++vertex) {
                const int index = bodyset.body[i].tri[k][vertex];
                if (index < 1 || index > bodyset.body[i].pointNum) {
                    fclose(infile);
                    throw std::runtime_error(
                        "Mesh vertex index out of range / 网格顶点编号越界: "
                        + file_utils::absolutePathForMessage(meshPath)
                        + ", face / 面 "
                        + std::to_string(k + 1));
                }
            }
        }
        fclose(infile);
        cout << "  Body " << i + 1 << " triangle count = " << bodyset.body[i].Size << endl;

        double dis;
        bodyset.body[i].sphereRadius = 0;
        for (int k = 0; k < bodyset.body[i].pointNum; ++k)
        {
            dis = bodyset.body[i].point[k].norm();
            bodyset.body[i].sphereRadius = bodyset.body[i].sphereRadius > dis ? bodyset.body[i].sphereRadius : dis;
        }
        cout << "  bodyset.body[i].sphereRadius = " << bodyset.body[i].sphereRadius << endl;
        vector3d a1, a2, center_tri,a3;
        bodyset.body[i].maxSize = 0;
        bodyset.body[i].Shape = new STLTetMesh[bodyset.body[i].Size]; // per-triangle mesh data
        for (int k = 0; k < bodyset.body[i].Size; ++k)
        {
            bodyset.body[i].Shape[k].v2 = bodyset.body[i].point[bodyset.body[i].tri[k][0]-1];
            bodyset.body[i].Shape[k].v3 = bodyset.body[i].point[bodyset.body[i].tri[k][1]-1];
            bodyset.body[i].Shape[k].v4 = bodyset.body[i].point[bodyset.body[i].tri[k][2]-1];
            
            a1 = bodyset.body[i].Shape[k].v3 - bodyset.body[i].Shape[k].v2;
            a2 = bodyset.body[i].Shape[k].v4 - bodyset.body[i].Shape[k].v2;
            a3 = bodyset.body[i].Shape[k].v3 - bodyset.body[i].Shape[k].v4;
            center_tri = (bodyset.body[i].Shape[k].v2 + bodyset.body[i].Shape[k].v3 + bodyset.body[i].Shape[k].v4) / 3;
            double meanSize = (a1.norm() + a2.norm() + a3.norm()) / 3;

            const vector3d areaNormal = a1.cross(a2); // outward normal, scaled by twice the triangle area
            bodyset.body[i].Shape[k].normal = areaNormal.normalized();
            bodyset.body[i].Shape[k].v1 = center_tri - bodyset.body[i].Shape[k].normal * meanSize;

            bodyset.body[i].maxSize = bodyset.body[i].maxSize > a1.norm() ? bodyset.body[i].maxSize : a1.norm();
            bodyset.body[i].maxSize = bodyset.body[i].maxSize > a2.norm() ? bodyset.body[i].maxSize : a2.norm();
            bodyset.body[i].maxSize = bodyset.body[i].maxSize > a3.norm() ? bodyset.body[i].maxSize : a3.norm();

            double judge = center_tri.dot(bodyset.body[i].Shape[k].normal);
            if (judge < 0) {
                std::cout << "Calculated outer normal error"<<std::endl;
            }
        } 
        std::cout << "  Max triangle size: " << bodyset.body[i].maxSize * 1000 << " mm" << std::endl;
    }

    for (int ib = 0; ib < bodyset.Num; ib++)	// allocate per-triangle AABB buffers
		bodyset.body[ib].eAABB = new AABB[bodyset.body[ib].Size];
	
	bodyset.SumFaceNum = 0;
    for (int ib = 0; ib < bodyset.Num; ib++) // total face count
		bodyset.SumFaceNum += bodyset.body[ib].Size;

    bodyset.IvDCM = new matrix3d[bodyset.Num];
    bodyset.DCM = new matrix3d[bodyset.Num];

    bodyset.triIdx = new int** [bodyset.Num];
    for (int i = 0; i < bodyset.Num; i++) {
        bodyset.triIdx[i] = new int* [bodyset.body[i].Size];
        for (int j = 0; j < bodyset.body[i].Size; j++)
        {
            bodyset.triIdx[i][j] = new int[6];
        }
    }
    bodyset.BodyLowerLmt = new vector3d[bodyset.Num];
    bodyset.BodyUpperLmt = new vector3d[bodyset.Num];

    std::cout << "LoadBodys Done......  \n" << endl;
}

void BODYSET::StateOutput(double time,double **force,int i,double **impulse)
{
	int n = Num;
	ofstream resfile;
	ostringstream convert;
	string Filename;

	Filename = "Data/DATA/OutputFile/state_bodys/At.";
	convert << setw(WidthOutput) << setfill('0') << i + 1;
	Filename.append(convert.str());
	Filename.append(".bt");
	//
	const auto outputPath = file_utils::prepareOutputFile(
		Filename, "rigid-body state output / 刚体状态输出");
	resfile.open(outputPath);
	if (!resfile.is_open()) {
		throw file_utils::pathError(
			"Cannot open / 无法打开",
			"rigid-body state output / 刚体状态输出",
			outputPath);
	}
	resfile << setiosflags(ios::scientific) << setprecision(PrecDouble); // output format
	//
    for (int s = 0; s < n; s++) // write body state
	{
			resfile << setw(WidthDouble) << time \
				<< setw(WidthDouble) << body[s].MassCenter[0] \
				<< setw(WidthDouble) << body[s].MassCenter[1] \
				<< setw(WidthDouble) << body[s].MassCenter[2] \
				<< setw(WidthDouble) << body[s].Vel[0] \
				<< setw(WidthDouble) << body[s].Vel[1] \
				<< setw(WidthDouble) << body[s].Vel[2] \
				<< setw(WidthDouble) << body[s].AngularVel[0] \
				<< setw(WidthDouble) << body[s].AngularVel[1] \
				<< setw(WidthDouble) << body[s].AngularVel[2] \
				<< setw(WidthDouble) << body[s].orien[0] \
				<< setw(WidthDouble) << body[s].orien[1] \
				<< setw(WidthDouble) << body[s].orien[2] \
				<< setw(WidthDouble) << body[s].orien[3] \
				<< setw(WidthDouble) << force[s][0]\
				<< setw(WidthDouble) << force[s][1]\
				<< setw(WidthDouble) << force[s][2]\
				<< setw(WidthDouble) << force[s][3]\
				<< setw(WidthDouble) << force[s][4]\
					<< setw(WidthDouble) << force[s][5];
			if (impulse) {
				for (int axis = 0; axis < 6; ++axis)
					resfile << setw(WidthDouble) << impulse[s][axis];
			}
			resfile << endl;
	}

	resfile.close();
	//
	convert.str("");
}


// Clean numeric cells (especially scientific notation with commas).
std::string cleanCellContent(const std::string& cell) {
    // Numeric cell: remove spaces and formatting characters.
    if (!cell.empty() && std::isdigit(static_cast<unsigned char>(cell[0]))) {
        std::string cleaned = cell;
        cleaned.erase(std::remove_if(cleaned.begin(), cleaned.end(),
            [](char c) { return std::isspace(c); }),
            cleaned.end());

        // Handle scientific notation.
        size_t ePos = cleaned.find('E');
        if (ePos != std::string::npos) {
            std::string mantissa = cleaned.substr(0, ePos);
            std::string exponent = cleaned.substr(ePos + 1);

            // Remove thousands separators.
            mantissa.erase(std::remove(mantissa.begin(), mantissa.end(), ','), mantissa.end());

            // Rebuild normalized string.
            cleaned = mantissa + "E" + exponent;
        }
        else {
            // Remove thousands separators.
            cleaned.erase(std::remove(cleaned.begin(), cleaned.end(), ','), cleaned.end());
        }

        return cleaned;
    }
    // Text cell: trim leading/trailing spaces only.
    else {
        // Trim leading spaces.
        size_t start = 0;
        while (start < cell.size() && std::isspace(static_cast<unsigned char>(cell[start]))) {
            ++start;
        }

        // Trim trailing spaces.
        size_t end = cell.size();
        while (end > start && std::isspace(static_cast<unsigned char>(cell[end - 1]))) {
            --end;
        }

        return cell.substr(start, end - start);
    }
}


// Split a CSV line into cleaned tokens.
inline std::vector<std::string> splitString(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::istringstream tokenStream(s);
    std::string token;

    while (std::getline(tokenStream, token, delimiter)) {
        token = cleanCellContent(token);
        if (!token.empty()) {
            tokens.push_back(token);
        }
    }

    return tokens;
}

// Load rigid bodies from CSV description.
inline bool LoadBodysFromCSV(BODYSET& bodyset, const std::string& filename) {
    const auto csvPath = file_utils::requireInputFile(
        filename, "rigid-body list file / 刚体列表文件");
    std::ifstream file(csvPath);
    if (!file.is_open()) {
        throw file_utils::pathError(
            "Cannot open / 无法打开",
            "rigid-body list file / 刚体列表文件",
            csvPath);
    }

    // Read first line (body count).
    std::string line;
    if (!std::getline(file, line)) {
        throw file_utils::pathError(
            "Invalid contents / 文件内容无效",
            "rigid-body list file / 刚体列表文件", csvPath,
            "file is empty or unreadable");
    }

    // Clean first line.
    line = cleanCellContent(line);

    int n;
    try {
        n = std::stoi(line);
        std::cout << "  Loaded " << n << " rigid bodies" << std::endl;
    }
    catch (const std::exception& e) {
        throw file_utils::pathError(
            "Invalid contents / 文件内容无效",
            "rigid-body list file / 刚体列表文件", csvPath,
            "invalid body count '" + line + "': " + e.what());
    }

    if (n <= 0) {
        throw file_utils::pathError(
            "Invalid contents / 文件内容无效",
            "rigid-body list file / 刚体列表文件", csvPath,
            "body count must be positive, got " + std::to_string(n));
    }

    // Read header line.
    if (!std::getline(file, line)) {
        throw file_utils::pathError(
            "Invalid contents / 文件内容无效",
            "rigid-body list file / 刚体列表文件", csvPath,
            "missing CSV header line");
    }

    // Allocate memory.
    bodyset.body = new BODY[n];
    bodyset.Num = n;

    // Read each body row.
    for (int i = 0; i < n; i++) {
        if (!std::getline(file, line)) {
            throw file_utils::pathError(
                "Invalid contents / 文件内容无效",
                "rigid-body list file / 刚体列表文件", csvPath,
                "missing row for body #" + std::to_string(i + 1));
        }

        // Split CSV row.
        auto tokens = splitString(line, ',');
        // Expect at least 24 columns (including name).
        if (tokens.size() < 24) {
            throw file_utils::pathError(
                "Invalid contents / 文件内容无效",
                "rigid-body list file / 刚体列表文件", csvPath,
                "body #" + std::to_string(i + 1) + " has "
                + std::to_string(tokens.size()) + " columns; at least 24 required");
        }

        try {
            // Read name (path).
            bodyset.body[i].name =
                file_utils::normalizePath(tokens[0]).u8string();

            // Read mass.
            bodyset.body[i].Mass = std::stod(tokens[1]);

            // Read mass center.
            bodyset.body[i].MassCenter[0] = std::stod(tokens[2]);
            bodyset.body[i].MassCenter[1] = std::stod(tokens[3]);
            bodyset.body[i].MassCenter[2] = std::stod(tokens[4]);

            // Read linear velocity.
            bodyset.body[i].Vel[0] = std::stod(tokens[5]);
            bodyset.body[i].Vel[1] = std::stod(tokens[6]);
            bodyset.body[i].Vel[2] = std::stod(tokens[7]);

            // Read angular velocity.
            bodyset.body[i].AngularVel[0] = std::stod(tokens[8]);
            bodyset.body[i].AngularVel[1] = std::stod(tokens[9]);
            bodyset.body[i].AngularVel[2] = std::stod(tokens[10]);

            // Read orientation quaternion.
            bodyset.body[i].orien[0] = std::stod(tokens[11]);
            bodyset.body[i].orien[1] = std::stod(tokens[12]);
            bodyset.body[i].orien[2] = std::stod(tokens[13]);
            bodyset.body[i].orien[3] = std::stod(tokens[14]);

            // Read inertia tensor.
            matrix3d inertia;
            for (int j = 0; j < 9; j++) {
                inertia[j / 3][j % 3] = std::stod(tokens[15 + j]);
            } 
            
            bodyset.body[i].I = inertia;

            // Optional trailing column keeps legacy 24-column files valid.
            if (tokens.size() >= 25) {
                std::size_t consumed = 0;
                bodyset.body[i].state = std::stoi(tokens[24], &consumed);
                if (consumed != tokens[24].size()) {
                    throw std::runtime_error("state is not an integer: " + tokens[24]);
                }
            }
            if (bodyset.body[i].state < BODY_DYNAMIC
                || bodyset.body[i].state > BODY_PRESCRIBED) {
                throw std::runtime_error(
                    "state must be 0 (dynamic), 1 (fixed), or 2 (prescribed)");
            }
            if (bodyset.body[i].state == BODY_FIXED) {
                bodyset.body[i].Vel = vector3d{0.0, 0.0, 0.0};
                bodyset.body[i].AngularVel = vector3d{0.0, 0.0, 0.0};
            }

            std::cout << "  Loaded body #" << i + 1 << ": " << bodyset.body[i].name << std::endl;
        }
        catch (const std::exception& e) {
            throw file_utils::pathError(
                "Invalid contents / 文件内容无效",
                "rigid-body list file / 刚体列表文件", csvPath,
                "cannot parse body #" + std::to_string(i + 1) + ": " + e.what());
        }
    }

    return true;
}
