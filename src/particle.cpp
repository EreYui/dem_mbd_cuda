#include "particle.h"
#include "force.h"
#include "constant.h"
#include "file_utils.h"
#include <sstream>
#include<iostream>
#include <fstream>  // std::ofstream
#include <iomanip>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include<vector>

using namespace std;

void PARTICLE::LoadParticles(string iniconfile)
{
    int j, count;
    string line;
    const auto particlePath = file_utils::requireInputFile(
        iniconfile, "particle input file / 颗粒输入文件");
    ifstream countfile(particlePath);
    if (!countfile.is_open()) {
        throw file_utils::pathError(
            "Cannot open / 无法打开",
            "particle input file / 颗粒输入文件",
            particlePath);
    }

    count = 0;
    while (getline(countfile, line))
    {
        if (line.find_first_not_of(" \t\r\n") != std::string::npos) {
            count++;
        }
    }
    countfile.close();

    if (count == 0) {
        throw file_utils::pathError(
            "Empty input / 输入文件为空",
            "particle input file / 颗粒输入文件",
            particlePath);
    }

    //Particle Alloc pt
    Num = count;
    Mass = new double[Num];
    Radius = new double[Num];
    Inertia = new vector3d[Num];
    Number = new int[Num];
    Status = new int[Num];

    Pos = new vector3d[Num];
    Vel = new vector3d[Num];
    Quat = new vector4d[Num];
    AngSpd = new vector3d[Num];

    Idx = new vector3i[Num];

    FILE* infile = file_utils::openInputFile(
        particlePath.u8string(), "particle input file / 颗粒输入文件");
    for (j = 0; j < Num; j++) {
        const int fieldsRead = fscanf(
            infile,
            "%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
            &Number[j], &Status[j], &Mass[j], &Radius[j],
            &Pos[j][0], &Pos[j][1], &Pos[j][2],
            &Vel[j][0], &Vel[j][1], &Vel[j][2],
            &Quat[j][0], &Quat[j][1], &Quat[j][2], &Quat[j][3],
            &AngSpd[j][0], &AngSpd[j][1], &AngSpd[j][2]);
        if (fieldsRead != 17) {
            fclose(infile);
            throw std::runtime_error(
                "Invalid particle record / 颗粒数据格式错误: "
                + file_utils::absolutePathForMessage(particlePath)
                + ", record / 记录 "
                + std::to_string(j + 1) + ", fields read / 已读取字段 "
                + std::to_string(fieldsRead) + "/17");
        }
    }
    fclose(infile);
    //
    for (j = 0; j < Num; j++)
    {
        const double sphereInertia = 0.4 * Mass[j] * Radius[j] * Radius[j];
        Inertia[j] = vector3d(sphereInertia, sphereInertia, sphereInertia);
    }

    MeshSize = 0.0;
    for (int i = 0; i < Num; i++)
        if (Radius[i] > MeshSize)	MeshSize = Radius[i];

    maxRadius = MeshSize;
    MeshSize = MeshSize * 2.5;

    cout << "  The number of particles  = " << Num << endl;
    cout << "  The radius of particles  = " << Radius[0] << endl;
    cout << "  The principal inertia of particles = "
         << Inertia[0][0] << ", " << Inertia[0][1] << ", " << Inertia[0][2]
         << endl;

    cout << "Particles Loaded......  \n" << endl;
}

void PARTICLE::LoadPrincipalInertia(const string& inertiafile)
{
    if (inertiafile.empty()) {
        return; // LoadParticles already synthesized the legacy solid-sphere value.
    }
    const auto path = file_utils::requireInputFile(
        inertiafile, "particle principal-inertia sidecar");
    ifstream input(path);
    if (!input) {
        throw file_utils::pathError(
            "Cannot open / 无法打开", "particle principal-inertia sidecar", path);
    }
    unordered_map<int, int> ownerRows;
    ownerRows.reserve(static_cast<size_t>(Num));
    for (int i = 0; i < Num; ++i) {
        ownerRows.emplace(Number[i], i);
    }
    vector<bool> seen(static_cast<size_t>(Num), false);
    string line;
    int lineNumber = 0;
    int records = 0;
    while (getline(input, line)) {
        ++lineNumber;
        if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
        istringstream record(line);
        int ownerId = 0;
        double ix = 0.0, iy = 0.0, iz = 0.0;
        string trailing;
        if (!(record >> ownerId >> ix >> iy >> iz) || (record >> trailing)) {
            throw runtime_error(
                "Invalid principal-inertia record at "
                + file_utils::absolutePathForMessage(path) + ":"
                + to_string(lineNumber) + "; expected owner_id Ixx Iyy Izz");
        }
        const auto owner = ownerRows.find(ownerId);
        if (owner == ownerRows.end()) {
            throw runtime_error(
                "Unknown owner ID in principal-inertia sidecar: "
                + to_string(ownerId));
        }
        const int row = owner->second;
        if (seen[static_cast<size_t>(row)]) {
            throw runtime_error(
                "Duplicate owner ID in principal-inertia sidecar: "
                + to_string(ownerId));
        }
        if (!isfinite(ix) || !isfinite(iy) || !isfinite(iz)
            || ix <= 0.0 || iy <= 0.0 || iz <= 0.0) {
            throw runtime_error(
                "Principal inertia must be finite and strictly positive for owner "
                + to_string(ownerId));
        }
        Inertia[row] = vector3d(ix, iy, iz);
        seen[static_cast<size_t>(row)] = true;
        ++records;
    }
    if (records != Num) {
        throw runtime_error(
            "Principal-inertia sidecar must contain exactly one record per owner");
    }
}

void PARTICLE::LoadClumpComponents(const string& componentfile)
{
    delete[] ComponentId;
    delete[] ComponentOwner;
    delete[] ComponentRadius;
    delete[] ComponentPosBody;
    ComponentId = nullptr;
    ComponentOwner = nullptr;
    ComponentRadius = nullptr;
    ComponentPosBody = nullptr;
    ComponentNum = 0;

    if (componentfile.empty()) {
        ComponentNum = Num;
        ComponentId = new int[ComponentNum];
        ComponentOwner = new int[ComponentNum];
        ComponentRadius = new double[ComponentNum];
        ComponentPosBody = new vector3d[ComponentNum];
        for (int i = 0; i < Num; ++i) {
            if (Number[i] < 0) {
                throw runtime_error(
                    "Default component ID requires a non-negative particle ID");
            }
            ComponentId[i] = Number[i];
            ComponentOwner[i] = i;
            ComponentRadius[i] = Radius[i];
            ComponentPosBody[i] = vector3d(0.0, 0.0, 0.0);
        }
        return;
    }

    const auto path = file_utils::requireInputFile(
        componentfile, "clump-component sidecar");
    ifstream countInput(path);
    string line;
    while (getline(countInput, line)) {
        if (line.find_first_not_of(" \t\r\n") != string::npos) ++ComponentNum;
    }
    if (ComponentNum <= 0) {
        throw runtime_error("Clump-component sidecar must not be empty");
    }
    ComponentId = new int[ComponentNum];
    ComponentOwner = new int[ComponentNum];
    ComponentRadius = new double[ComponentNum];
    ComponentPosBody = new vector3d[ComponentNum];

    unordered_map<int, int> ownerRows;
    ownerRows.reserve(static_cast<size_t>(Num));
    for (int i = 0; i < Num; ++i) ownerRows.emplace(Number[i], i);
    unordered_set<int> componentIds;
    componentIds.reserve(static_cast<size_t>(ComponentNum));
    vector<int> ownerCounts(static_cast<size_t>(Num), 0);
    ifstream input(path);
    int row = 0;
    int lineNumber = 0;
    while (getline(input, line)) {
        ++lineNumber;
        if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
        istringstream record(line);
        int componentId = 0, ownerId = 0;
        double radius = 0.0, x = 0.0, y = 0.0, z = 0.0;
        string trailing;
        if (!(record >> componentId >> ownerId >> radius >> x >> y >> z)
            || (record >> trailing)) {
            throw runtime_error(
                "Invalid clump-component record at "
                + file_utils::absolutePathForMessage(path) + ":"
                + to_string(lineNumber)
                + "; expected component_id owner_id radius body_x body_y body_z");
        }
        const auto owner = ownerRows.find(ownerId);
        if (owner == ownerRows.end()) {
            throw runtime_error(
                "Unknown owner ID in clump-component sidecar: "
                + to_string(ownerId));
        }
        if (!componentIds.insert(componentId).second) {
            throw runtime_error(
                "Duplicate component ID in clump-component sidecar: "
                + to_string(componentId));
        }
        if (componentId < 0) {
            throw runtime_error("Clump component IDs must lie in [0, INT32_MAX]");
        }
        if (!isfinite(radius) || radius <= 0.0
            || !isfinite(x) || !isfinite(y) || !isfinite(z)) {
            throw runtime_error(
                "Clump-component geometry must be finite with positive radius");
        }
        ComponentId[row] = componentId;
        ComponentOwner[row] = owner->second;
        ComponentRadius[row] = radius;
        ComponentPosBody[row] = vector3d(x, y, z);
        ++ownerCounts[static_cast<size_t>(owner->second)];
        ++row;
    }
    for (int i = 0; i < Num; ++i) {
        if (ownerCounts[static_cast<size_t>(i)] == 0) {
            throw runtime_error(
                "Every particle owner must have at least one clump component; missing owner "
                + to_string(Number[i]));
        }
        if (ownerCounts[static_cast<size_t>(i)] > 8) {
            throw runtime_error(
                "Each particle owner may have at most 8 clump components; owner "
                + to_string(Number[i]));
        }
    }
}

void PARTICLE::FreeParticles()
{
    delete[]Mass;
    delete[]Radius;
    delete[]Inertia;
    delete[]Number;
    delete[]Status;
    delete[]Pos;
    delete[]Vel;
    delete[]Quat;
    delete[]AngSpd;
    delete[]Idx;
    delete[]ComponentId;
    delete[]ComponentOwner;
    delete[]ComponentRadius;
    delete[]ComponentPosBody;
    Mass = nullptr;
    Radius = nullptr;
    Inertia = nullptr;
    Number = nullptr;
    Status = nullptr;
    Pos = nullptr;
    Vel = nullptr;
    Quat = nullptr;
    AngSpd = nullptr;
    Idx = nullptr;
    ComponentId = nullptr;
    ComponentOwner = nullptr;
    ComponentRadius = nullptr;
    ComponentPosBody = nullptr;
    ComponentNum = 0;
    Num = 0;
}

PARTICLE::~PARTICLE()
{
    FreeParticles();
}


void PARTICLE::StateOutput(int i)
{
	ofstream resfile;
	ostringstream convert;
	string Filename;

	Filename = "Data/DATA/OutputFile/state_particles/ph.";
	convert << setw(WidthOutput) << setfill('0') << i + 1;
	Filename.append(convert.str());
	Filename.append(".bt");
	//
	const int WidthInt = 10;
	const int WidthDouble = 16;
	const int PrecDouble = 12;

	int num_threads = 1;
	std::vector<std::string> buffers(num_threads);

	{
		int tid = 0;
		// 预分配每线程缓冲区
		buffers[tid].reserve((WidthInt * 2 + WidthDouble * 15 + 1) * (Num / num_threads));
	}

	for (int j = 0; j < Num; ++j) {
		int tid = 0;
		char line[1024];
		snprintf(line, sizeof(line),
			"%*d %*d %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e\n",
			WidthInt, Number[j],
			WidthInt, Status[j],
			WidthDouble, PrecDouble, Mass[j],
			WidthDouble, PrecDouble, Radius[j],
			WidthDouble, PrecDouble, Pos[j][0],
			WidthDouble, PrecDouble, Pos[j][1],
			WidthDouble, PrecDouble, Pos[j][2],
			WidthDouble, PrecDouble, Vel[j][0],
			WidthDouble, PrecDouble, Vel[j][1],
			WidthDouble, PrecDouble, Vel[j][2],
			WidthDouble, PrecDouble, Quat[j][0],
			WidthDouble, PrecDouble, Quat[j][1],
			WidthDouble, PrecDouble, Quat[j][2],
			WidthDouble, PrecDouble, Quat[j][3],
			WidthDouble, PrecDouble, AngSpd[j][0],
			WidthDouble, PrecDouble, AngSpd[j][1],
			WidthDouble, PrecDouble, AngSpd[j][2]
		);
		buffers[tid] += line;
	}

	const auto outputPath = file_utils::prepareOutputFile(
		Filename, "particle state output / 颗粒状态输出");
	std::ofstream file(outputPath);
	if (!file.is_open()) {
		throw file_utils::pathError(
			"Cannot open / 无法打开",
			"particle state output / 颗粒状态输出",
			outputPath);
	}

	for (auto& buf : buffers) {
		file << buf;
	}

}

//输出颗粒受到的力,算例中包含刚体
void FORCE::ParticlesForceOutput(int i, int nb, int nbody)
{
	ofstream resfile;
	ostringstream convert;
	string Filename;

	Filename = "Data/DATA/OutputFile/force_particles/F.";
	convert << setw(WidthOutput) << setfill('0') << i + 1;
	Filename.append(convert.str());
	Filename.append(".bt");
	//
	const int WidthInt = 10;
	const int WidthDouble = 16;
	const int PrecDouble = 12;

	int num_threads = 1;
	std::vector<std::string> buffers(num_threads);

	{
		int tid = 0;
		// 预分配每线程缓冲区
		buffers[tid].reserve((WidthDouble * (21 + nbody * 6)) * (nb / num_threads));
	}

	for (int j = 0; j < nb; ++j) {
		int tid = 0;
		char line[32768];// 增大缓冲区防止溢出
		int offset = 0;
		// 生成前21个数据
		//Fc,颗粒间的接触力,6
		//Fw,墙体给颗粒的力，6
		//Ff,万有引力和重力，3
		//FR,颗粒受到的所有的力的合力，6
		offset += snprintf(line + offset, sizeof(line) - offset,
			"%*.*e %*.*e %*.*e %*.*e %*.*e %*.*e "
			"%*.*e %*.*e %*.*e %*.*e %*.*e %*.*e "
			"%*.*e %*.*e %*.*e "
			"%*.*e %*.*e %*.*e %*.*e %*.*e %*.*e",
			WidthDouble, PrecDouble, Fc[j][0],
			WidthDouble, PrecDouble, Fc[j][1],
			WidthDouble, PrecDouble, Fc[j][2],
			WidthDouble, PrecDouble, Fc[j][3],
			WidthDouble, PrecDouble, Fc[j][4],
			WidthDouble, PrecDouble, Fc[j][5],
			WidthDouble, PrecDouble, Fw[j][0],
			WidthDouble, PrecDouble, Fw[j][1],
			WidthDouble, PrecDouble, Fw[j][2],
			WidthDouble, PrecDouble, Fw[j][3],
			WidthDouble, PrecDouble, Fw[j][4],
			WidthDouble, PrecDouble, Fw[j][5],
			WidthDouble, PrecDouble, Ff[j][0],
			WidthDouble, PrecDouble, Ff[j][1],
			WidthDouble, PrecDouble, Ff[j][2],
			WidthDouble, PrecDouble, FR[j][0],
			WidthDouble, PrecDouble, FR[j][1],
			WidthDouble, PrecDouble, FR[j][2],
			WidthDouble, PrecDouble, FR[j][3],
			WidthDouble, PrecDouble, FR[j][4],
			WidthDouble, PrecDouble, FR[j][5]
		);
		// append nbody*6 FI_pt values
		for (int k = 0; k < nbody; ++k) {
			for (int m = 0; m < 6; ++m) {
				offset += snprintf(line + offset, sizeof(line) - offset,
					" %*.*e", WidthDouble, PrecDouble, FI_pt[k][j][m]);
			}
		}
		// append newline
		offset += snprintf(line + offset, sizeof(line) - offset, "\n");

		buffers[tid] += line;
	}

	const auto outputPath = file_utils::prepareOutputFile(
		Filename, "particle force output / 颗粒受力输出");
	std::ofstream file(outputPath);
	if (!file.is_open()) {
		throw file_utils::pathError(
			"Cannot open / 无法打开",
			"particle force output / 颗粒受力输出",
			outputPath);
	}

	for (auto& buf : buffers) {
		file << buf;
	}
}

//输出颗粒受到的力,算例中无刚体
void FORCE::ParticlesForceOutput(int i, int nb)
{
	ofstream resfile;
	ostringstream convert;
	string Filename;

	Filename = "Data/DATA/OutputFile/force_particles/F.";
	convert << setw(WidthOutput) << setfill('0') << i + 1;
	Filename.append(convert.str());
	Filename.append(".bt");
	//
	const int WidthInt = 10;
	const int WidthDouble = 16;
	const int PrecDouble = 12;

	int num_threads = 1;
	std::vector<std::string> buffers(num_threads);

	{
		int tid = 0;
		// 预分配每线程缓冲区
		buffers[tid].reserve((WidthDouble * 21 + 21) * (nb / num_threads));
	}

	for (int j = 0; j < nb; ++j) {
		int tid = 0;
		char line[1024];
		snprintf(line, sizeof(line),
			"%*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e %*.*e\n",
			WidthDouble, PrecDouble, Fc[j][0],
			WidthDouble, PrecDouble, Fc[j][1],
			WidthDouble, PrecDouble, Fc[j][2],
			WidthDouble, PrecDouble, Fc[j][3],
			WidthDouble, PrecDouble, Fc[j][4],
			WidthDouble, PrecDouble, Fc[j][5],
			WidthDouble, PrecDouble, Fw[j][0],
			WidthDouble, PrecDouble, Fw[j][1],
			WidthDouble, PrecDouble, Fw[j][2],
			WidthDouble, PrecDouble, Fw[j][3],
			WidthDouble, PrecDouble, Fw[j][4],
			WidthDouble, PrecDouble, Fw[j][5],
			WidthDouble, PrecDouble, Ff[j][0],
			WidthDouble, PrecDouble, Ff[j][1],
			WidthDouble, PrecDouble, Ff[j][2],
			WidthDouble, PrecDouble, FR[j][0],
			WidthDouble, PrecDouble, FR[j][1],
			WidthDouble, PrecDouble, FR[j][2],
			WidthDouble, PrecDouble, FR[j][3],
			WidthDouble, PrecDouble, FR[j][4],
			WidthDouble, PrecDouble, FR[j][5]
		);
		buffers[tid] += line;
	}

	const auto outputPath = file_utils::prepareOutputFile(
		Filename, "particle force output / 颗粒受力输出");
	std::ofstream file(outputPath);
	if (!file.is_open()) {
		throw file_utils::pathError(
			"Cannot open / 无法打开",
			"particle force output / 颗粒受力输出",
			outputPath);
	}

	for (auto& buf : buffers) {
		file << buf;
	}
}

//输出刚体受到的力
void FORCE::BodysetForceOutput(int i, int nb, int nbody)
{
	const int WidthOutput = 5;      // File index width
	const int WidthDouble = 16;     // Double output width
	const int PrecDouble = 12;     // Double output precision

	// Output file path
	std::ostringstream filenameStream;
	filenameStream << "Data/DATA/OutputFile/force_bodys/FB."
		<< std::setw(WidthOutput) << std::setfill('0') << (i + 1)
		<< ".bt";
	std::string Filename = filenameStream.str();

	// Get the number of threads
	int num_threads = 1;
	std::vector<std::ostringstream> threadBuffers(num_threads);

	// Total lines to output: number of bodies * number of grains per body
	const int totalLines = nbody * nb;

	{
		int tid = 0;
		// Pre-allocate buffer size for each thread, approximately 300 characters per line
		threadBuffers[tid].str().reserve(totalLines / num_threads * 350);
	}

	for (int index = 0; index < totalLines; ++index) {
		int tid = 0;
		int body_idx = index / nb;   // Body index
		int grain_idx = index % nb;   // Grain index

		// Get current body-grain forces
		//double* forces = ForceBody[body_idx][grain_idx];

		// Format output string
		for (int k = 0; k < 19; ++k) {
			threadBuffers[tid] << std::setw(WidthDouble)
				<< std::setprecision(PrecDouble)
				<< std::scientific
				<< ForceBody[body_idx][grain_idx][k];

			// Add space between values except after the last one
			if (k < 18) threadBuffers[tid] << " ";
		}
		threadBuffers[tid] << "\n"; // append newline
	}

	// Write to file
	const auto outputPath = file_utils::prepareOutputFile(
		Filename, "body force output / 刚体受力输出");
	std::ofstream outFile(outputPath, std::ios::binary);
	if (!outFile.is_open()) {
		throw file_utils::pathError(
			"Cannot open / 无法打开",
			"body force output / 刚体受力输出",
			outputPath);
	}

	// Merge thread buffers and write to file
	for (auto& buf : threadBuffers) {
		outFile << buf.str();
	}
}
