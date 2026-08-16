#include <algorithm>
#include <cassert>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "simulation.h"
#include "particle.h"
#include "multibody.h"
#include "grid.h"
#include "force.h"
#include "wall.h"
#include "state.h"
#include "dyneq.h"
#include "file_utils.h"

using namespace std;

namespace {

bool assignMechanicalParameter(
    const string& name, const string& prefix, double value, MECH& mech)
{
    if (name.compare(0, prefix.size(), prefix) != 0) {
        return false;
    }

    const string suffix = name.substr(prefix.size());
    if (suffix == "epsS") mech.epsS = value;
    else if (suffix == "mu") mech.mu = value;
    else if (suffix == "mu_T") mech.mu_T = value;
    else if (suffix == "mu_R") mech.mu_R = value;
    else if (suffix == "epsN") mech.epsN = value;
    else if (suffix == "cohesive") mech.c = value;
    else if (suffix == "beta") mech.beta = value;
    else if (suffix == "kN") mech.kN = value;
    else return false;
    return true;
}

bool assignLegacyMechanicalParameter(
    const string& name, double value, MECH& mechPP, MECH& mechPT, MECH& mechPW)
{
    if (!assignMechanicalParameter(name, "MECH_", value, mechPP)) {
        return false;
    }

    // A legacy MECH_* parameter applies to all contact types.
    assignMechanicalParameter(name, "MECH_", value, mechPT);
    assignMechanicalParameter(name, "MECH_", value, mechPW);
    return true;
}

bool parseOutputSwitch(
    const string& name, const string& value, const string& parameterFile, int lineNumber)
{
    if (value == "0") return false;
    if (value == "1") return true;
    throw runtime_error(
        name + " must be 0 or 1 / 必须为 0 或 1 at "
        + parameterFile + ":" + to_string(lineNumber));
}

void printMechanicalSet(const char* label, const MECH& mech)
{
    cout << "  " << label << ".epsS = " << mech.epsS << endl;
    cout << "  " << label << ".mu   = " << mech.mu << endl;
    cout << "  " << label << ".mu_T = " << mech.mu_T << endl;
    cout << "  " << label << ".mu_R = " << mech.mu_R << endl;
    cout << "  " << label << ".epsN = " << mech.epsN << endl;
    cout << "  " << label << ".c    = " << mech.c << endl;
    cout << "  " << label << ".beta = " << mech.beta << endl;
    cout << "  " << label << ".kN   = " << mech.kN << endl;
}

void writeMechanicalSet(ostream& out, const char* prefix, const MECH& mech)
{
    out << prefix << "epsS     = " << left << setw(8) << mech.epsS << endl;
    out << prefix << "mu       = " << left << setw(8) << mech.mu << endl;
    out << prefix << "mu_T     = " << left << setw(8) << mech.mu_T << endl;
    out << prefix << "mu_R     = " << left << setw(8) << mech.mu_R << endl;
    out << prefix << "epsN     = " << left << setw(8) << mech.epsN << endl;
    out << prefix << "cohesive = " << left << setw(8) << mech.c << endl;
    out << prefix << "beta     = " << left << setw(8) << mech.beta << endl;
    out << prefix << "kN       = " << left << setw(8) << mech.kN << endl;
}

} // namespace

Simulation::Simulation() : time(0) {
    // 初始化成员变量
}
Simulation::~Simulation() = default;


void Simulation::init(const std::string& parafile) {
    // 加载配置文件
    loadParas(parafile);

	// 加载颗粒数据
    if (control.ParticlesFile.empty()) {
        throw std::runtime_error(
            "CONTROL_PARTICLES_FILE is missing / 未配置颗粒文件路径");
    }
    pt.LoadParticles(control.ParticlesFile);
    
	// 加载刚体数据
    if (control.Multibody_flag) {
        if (control.MultibodyFile.empty()) {
            throw std::runtime_error(
                "CONTROL_MULTIBODY_FILE is missing while multibody is enabled"
                " / 已启用多刚体但未配置刚体文件路径");
        }
        bodyset.LoadBodys(control.MultibodyFile);
	}
	// 加载壁面数据
    if (control.Wall_flag) {
        if (control.WallFile.empty()) {
            throw std::runtime_error(
                "CONTROL_WALL_FILE is missing while walls are enabled"
                " / 已启用墙体但未配置墙体文件路径");
        }
        control.wall.LoadWalls(control.WallFile);
    }
}


void Simulation::run() {

    cout << "Start integrate......  " << endl;
    clock_t systim1;

    systim1 = clock(); 
    if (control.Multibody_flag)
    {
        IntegrateDemMultiBody();
        systim1 = clock() - systim1;
    }
    else
    {
        IntegrateDem();
        systim1 = clock() - systim1;
    }
    cout << "Time = " << systim1/60000<<" min" << endl;
    cout << "Integrate Done......  " << endl;

    int total_time = systim1 / 60000;
    if (control.OutputSummary) {
        summaryOutput(total_time);
    }
}

void Simulation::loadParas(const std::string& parafile)
{
    int numchar, id, n;
    int flag;

    ifstream infile;
    string line, core, name, value;

    flag = 0;

    const auto parameterPath = file_utils::requireInputFile(
        parafile, "simulation parameter file / 仿真参数文件");
    infile.open(parameterPath);
    if (!infile.is_open()) {
        throw file_utils::pathError(
            "Cannot open / 无法打开",
            "simulation parameter file / 仿真参数文件",
            parameterPath);
    }

    int lineNumber = 0;
    while (getline(infile, line))
    {
        ++lineNumber;
        numchar = int(line.length());
        for (n = 0; n < numchar; n++) if (line[n] == '#') break;
        core.assign(line, 0, n);
        //
        core.erase(std::remove(core.begin(), core.end(), '\t'), core.end());
        core.erase(std::remove(core.begin(), core.end(), ' '), core.end());
        //
        if (!core.empty())
        {
            //
            numchar = int(core.length());
            id = int(core.find('='));
            if (!(id > 0 && id + 1 < numchar)) {
                throw std::runtime_error(
                    "Invalid parameter line / 参数行格式错误: "
                    + file_utils::absolutePathForMessage(parameterPath) + ":"
                    + std::to_string(lineNumber));
            }
            name.assign(core, 0, id);
            n = numchar - id - 1;
            value.assign(core, id + 1, n);
            const double numericValue = atof(value.c_str());
            if (name == "ODE_Option") { 
                ode.Option = numericValue;                  flag++;
            } else if (name == "ODE_StepSize") {
                ode.StepSize = numericValue;                flag++;
            } else if (name == "ODE_EndStep") { 
                ode.EndStep = numericValue;                 flag++;
            } else if (name == "ODE_StartStep") { 
                ode.StartStep = numericValue;               flag++;
            } else if (name == "ODE_OutputInterval") { 
                ode.OutputInterval = numericValue;          flag++;
            } else if (
                assignMechanicalParameter(name, "MECH_PP_", numericValue, mechPP)
                || assignMechanicalParameter(name, "MECH_PT_", numericValue, mechPT)
                || assignMechanicalParameter(name, "MECH_PW_", numericValue, mechPW)
                || assignLegacyMechanicalParameter(name, numericValue, mechPP, mechPT, mechPW)) {
                flag++;
            } else if (name == "CONTROL_GRAVITY_FLAG") { 
                control.Grav_flag = numericValue;           flag++;
            } else if (name == "CONTROL_GRAVITY_X") { 
                control.g[0] = atof(value.c_str());         flag++;
            } else if (name == "CONTROL_GRAVITY_Y") { 
                control.g[1] = atof(value.c_str());         flag++;
            } else if (name == "CONTROL_GRAVITY_Z") { 
                control.g[2] = atof(value.c_str());         flag++;
            } else if (name == "CONTROL_UNIVERSAL_GRAVITATION_FLAG") { 
                control.Universal_gravitation_flag = atof(value.c_str());    flag++; 
            } else if (name == "CONTROL_WALL_FLAG") { 
                control.Wall_flag = atof(value.c_str());    flag++; 
            } else if (name == "CONTROL_WALL_FILE") { 
                control.WallFile = value;                   flag++; 
            }else if (name == "CONTROL_MULTIBODY_FLAG") { 
                control.Multibody_flag = atof(value.c_str());  flag++; 
            } else if (name == "CONTROL_MULTIBODY_FILE") { 
                control.MultibodyFile = value;                 flag++; 
            } else if (name == "CONTROL_PARTICLES_FILE") { 
                control.ParticlesFile = value;                 flag++; 
            } else if (name == "CONTROL_OUTPUT_PARTICLE_STATE") {
                control.OutputParticleState = parseOutputSwitch(
                    name, value, parameterPath.u8string(), lineNumber); flag++;
            } else if (name == "CONTROL_OUTPUT_PARTICLE_FORCE") {
                control.OutputParticleForce = parseOutputSwitch(
                    name, value, parameterPath.u8string(), lineNumber); flag++;
            } else if (name == "CONTROL_OUTPUT_BODY_STATE") {
                control.OutputBodyState = parseOutputSwitch(
                    name, value, parameterPath.u8string(), lineNumber); flag++;
            } else if (name == "CONTROL_OUTPUT_BODY_FORCE") {
                control.OutputBodyForce = parseOutputSwitch(
                    name, value, parameterPath.u8string(), lineNumber); flag++;
            } else if (name == "CONTROL_OUTPUT_SUMMARY") {
                control.OutputSummary = parseOutputSwitch(
                    name, value, parameterPath.u8string(), lineNumber); flag++;
            } else { 
                throw std::runtime_error(
                    "Unrecognized parameter / 无法识别的参数 '" + name
                    + "' at " + parameterPath.u8string() + ":"
                    + std::to_string(lineNumber));
            }
        }
    }
    infile.close();

    cout << "Loaded Parameters: " << endl;
    cout << "  ode.Option         = " << ode.Option << endl;
    cout << "  ode.StepSize       = " << ode.StepSize << endl;
    cout << "  ode.EndStep        = " << ode.EndStep << endl;
    cout << "  ode.StartStep      = " << ode.StartStep << endl;
    cout << "  ode.OutputInterval = " << ode.OutputInterval << endl;
    printMechanicalSet("mechPP", mechPP);
    printMechanicalSet("mechPT", mechPT);
    printMechanicalSet("mechPW", mechPW);
    cout << "  control.Grav_flag  = " << control.Grav_flag << endl;
    cout << "  control.g[0]       = " << control.g[0] << endl;
    cout << "  control.g[1]       = " << control.g[1] << endl;
    cout << "  control.g[2]       = " << control.g[2] << endl;
    cout << "  control.Wall_flag  = " << control.Wall_flag << endl;
    cout << "  control.WallFile   = " << control.WallFile << endl;
    cout << "  control.Multibody_flag = " << control.Multibody_flag << endl;
    cout << "  control.MultibodyFile  = " << control.MultibodyFile << endl;
    cout << "  control.ParticlesFile  = " << control.ParticlesFile << endl;
    cout << "  control.OutputParticleState = " << control.OutputParticleState << endl;
    cout << "  control.OutputParticleForce = " << control.OutputParticleForce << endl;
    cout << "  control.OutputBodyState     = " << control.OutputBodyState << endl;
    cout << "  control.OutputBodyForce     = " << control.OutputBodyForce << endl;
    cout << "  control.OutputSummary       = " << control.OutputSummary << endl;
    //cout<<"dahsdf"<<endl;
  /*  if (flag != NumParameter)
    {
        cout<<flag<<endl;
        cout<<"Important parameter setting missing!"<<endl;
        exit(1);
    }*/
    //
    cout << "Parameters Loaded......  \n" << endl;
}

#ifndef DEM_MBD_USE_CUDA
void Simulation::IntegrateDem()
{
    int count = 0;
    double C = 0.08 * pt.MeshSize;
    double dt = ode.StepSize / 2;
    FORCE force(pt.Num, 0);
    pt.DRmax = 0;

    //defined the contact pair:Cp  Contact pair
    Cp** part;
    Cp** poly_part;
    part = new Cp * [pt.Num];
    for (int i = 0; i < pt.Num; i++)
    {
        part[i] = (Cp*)malloc(sizeof(Cp));
        for (int k = 0; k < 3; k++)
            for (int n = 0; n < 3; n++)
                part[i]->deltP[k][n] = 0;
        part[i]->flag = 0;
        part[i]->idx = -1;
        part[i]->next = nullptr;
    }

    poly_part = new Cp * [pt.Num]; 
    for (int i = 0; i < pt.Num; i++)
    {
        poly_part[i] = (Cp*)malloc(sizeof(Cp));
        for (int k = 0; k < 3; k++)
            for (int n = 0; n < 3; n++)
                poly_part[i]->deltP[k][n] = 0;
        poly_part[i]->flag = 0;
        poly_part[i]->idx = -1;
        poly_part[i]->body = -1;
        poly_part[i]->next = nullptr;
    }

    VHALF vhalf_pt(pt.Num);
    ACC ptacc(pt.Num);

    CreatTree(pt, grid);
    force.ForceCaculate(0.0, pt, control, mechPP, mechPW, ode, part, poly_part, grid);
    if (control.OutputParticleForce) {
        force.ParticlesForceOutput(-1, pt.Num); // output initial force
    }
    if (control.OutputParticleState) {
        pt.StateOutput(-1); // output initial particle state
    }

    for (int j = 0; j < pt.Num; j++)
    {
        DynEqnParticle(j, pt, force.FR[j], ptacc.acc[j], ptacc.domg[j]);
        vhalf_pt.Vel[j] = pt.Vel[j] + ptacc.acc[j] * dt ;
        vhalf_pt.AngSpd[j] = pt.AngSpd[j] + ptacc.domg[j] * dt;
        vhalf_pt.dQuat[j] = gama(pt.Quat[j], vhalf_pt.AngSpd[j]);
    }

    int th_num = 1;//
    double* Vmax_thread;
    Vmax_thread = new double[th_num];
//---------------------------------------------------------------------------------------------------------------------
    cout << "Calculating......." << endl;
    for (int i = ode.StartStep; i < ode.EndStep; i++)
    {
        for (int j = 0; j < th_num; j++)Vmax_thread[j] = 0;
        double time = i * ode.StepSize;

        for (int j = 0; j < pt.Num; j++)
        {
            int cur_th_num = 0;//
            if (pt.Status[j] != -1)//integral
            {
                pt.Pos[j] += vhalf_pt.Vel[j] * dt * 2;
                pt.Quat[j] += vhalf_pt.dQuat[j] * dt * 2;
                pt.Quat[j] = pt.Quat[j].normalize();

                pt.Vel[j] = vhalf_pt.Vel[j] + ptacc.acc[j] * dt;
                pt.AngSpd[j] = vhalf_pt.AngSpd[j] + ptacc.domg[j] * dt;

                double P_tmp = (vhalf_pt.Vel[j] * dt * 2).norm();
                if (P_tmp > Vmax_thread[cur_th_num])Vmax_thread[cur_th_num] = P_tmp;//
            }
        }
        pt.Vmax = 0.0;//

        for (int j = 0; j < th_num; j++)
        {
            if (Vmax_thread[j] > pt.Vmax)
                pt.Vmax = Vmax_thread[j];
        }
        pt.DRmax = pt.DRmax + pt.Vmax;//tree

        if (pt.DRmax > C)//update mesh
        {
            CreatTree(pt, grid);
            pt.DRmax = 0;
        }

        force.ForceCaculate(time, pt, control, mechPP, mechPW, ode, part, poly_part, grid);
        for (int j = 0; j < pt.Num; j++)
        {
            DynEqnParticle(j, pt, force.FR[j], ptacc.acc[j], ptacc.domg[j]);
            vhalf_pt.Vel[j] += ptacc.acc[j] * 2 * dt;
            vhalf_pt.AngSpd[j] += ptacc.domg[j] * 2 * dt;
            vhalf_pt.dQuat[j] = gama(pt.Quat[j], vhalf_pt.AngSpd[j]);
        }

        // output check
        count++;
        if (ode.OutputInterval == count || i == ode.EndStep-1)//CT.FoutInterval
        {
            if (control.OutputParticleForce) {
                force.ParticlesForceOutput(i, pt.Num); // 输出颗粒受力
            }
            if (control.OutputParticleState) {
                pt.StateOutput(i); // 输出颗粒状态
            }
            count = 0;
            cout << "Step = " << i + 1 << "/" << ode.EndStep << ", Size: " << grid.size() << endl;
        }
    }//ode.EndStep

    delete[] part;
    delete[] poly_part;
    delete[] Vmax_thread;
}

// Integrate Dem Multi Body System
void Simulation::IntegrateDemMultiBody()
{
    bodyset.updateMatrix();

    int count = 0;
    double C = 0.08 * pt.MeshSize;
    FORCE force(pt.Num, bodyset.Num);
    double dt = ode.StepSize / 2;    
    pt.DRmax = 0;

    //defined the contact pair:Cp  Contact pair
    Cp** part;
    Cp** poly_part;
    part = new Cp * [pt.Num];
    for (int i = 0; i < pt.Num; i++)
    {
        part[i] = (Cp*)malloc(sizeof(Cp));
        for (int k = 0; k < 3; k++)
            for (int n = 0; n < 3; n++)
                part[i]->deltP[k][n] = 0;
        part[i]->flag = 0;
        part[i]->idx = -1;
        part[i]->body = -1;
        part[i]->next = nullptr;
    }

    poly_part = new Cp * [pt.Num]; 
    for (int i = 0; i < pt.Num; i++)
    {
        poly_part[i] = (Cp*)malloc(sizeof(Cp));
        for (int k = 0; k < 3; k++)
            for (int n = 0; n < 3; n++)
                poly_part[i]->deltP[k][n] = 0;
        poly_part[i]->flag = 0;
        poly_part[i]->idx = -1;
        poly_part[i]->body = -1;
        poly_part[i]->next = nullptr;
    }

    VHALF vhalf_bd(bodyset.Num), vhalf_pt(pt.Num);
    ACC bodyacc(bodyset.Num), ptacc(pt.Num);

    CreatTree(pt, bodyset, grid);
    force.ForceCaculate(
        0.0, pt, bodyset, control, mechPP, mechPT, mechPW, ode,
        part, poly_part, grid, vhalf_pt, vhalf_bd);

    if (control.OutputParticleForce) {
        force.ParticlesForceOutput(-1, pt.Num, bodyset.Num);
    }
    if (control.OutputBodyForce) {
        force.BodysetForceOutput(-1, pt.Num, bodyset.Num);
    }
    if (control.OutputParticleState) {
        pt.StateOutput(-1);
    }
    if (control.OutputBodyState) {
        bodyset.StateOutput(0.0, force.poly_F, -1);
    }

    double** ddu = new double* [bodyset.Num];
    for (int i = 0; i < bodyset.Num; i++) {
        ddu[i] = new double[6]();
    }
    for (int j = 0; j < pt.Num; j++)
    {
        DynEqnParticle(j, pt, force.FR[j], ptacc.acc[j], ptacc.domg[j]);
        vhalf_pt.Vel[j] = pt.Vel[j] + ptacc.acc[j] * dt;
        vhalf_pt.AngSpd[j] = pt.AngSpd[j] + ptacc.domg[j] * dt;
        vhalf_pt.quat[j] = (pt.Quat[j] + gama(pt.Quat[j], pt.AngSpd[j]) * dt).normalize();
        vhalf_pt.dQuat[j] = gama(vhalf_pt.quat[j], vhalf_pt.AngSpd[j]);
    }

    DynEqnBodySet(0.0, bodyset, force.poly_F, ddu);
    for (int j = 0; j < bodyset.Num; j++)
    {
        bodyacc.acc[j][0] = ddu[j][0]; bodyacc.acc[j][1] = ddu[j][1]; bodyacc.acc[j][2] = ddu[j][2];
        bodyacc.domg[j][0] = ddu[j][3]; bodyacc.domg[j][1] = ddu[j][4]; bodyacc.domg[j][2] = ddu[j][5];

        vhalf_bd.Vel[j] = bodyset.body[j].Vel + bodyacc.acc[j] * dt;
        vhalf_bd.AngSpd[j] = bodyset.body[j].AngularVel + bodyacc.domg[j] * dt;
        vhalf_bd.quat[j] = (bodyset.body[j].orien + gama(bodyset.body[j].orien, bodyset.body[j].AngularVel) * dt).normalize();
        vhalf_bd.dQuat[j] = gama(vhalf_bd.quat[j], vhalf_bd.AngSpd[j]);
    }
    //-------------------------------------------------------------------------------------------------------------------------------------
    int th_num = 1;//
    double* Vmax_thread;
    Vmax_thread = new double[th_num];

    cout << "Calculating......." << endl;
    for (int i = ode.StartStep; i < ode.EndStep; i++)
    {
        for (int j = 0; j < th_num; j++)    Vmax_thread[j] = 0.0;
        double time = i * ode.StepSize;

        for (int j = 0; j < pt.Num; j++)
        {
            int cur_th_num = 0;//
            if (pt.Status[j] != -1)//integral
            {
                pt.Pos[j] += vhalf_pt.Vel[j] * dt * 2;
                pt.Quat[j] += vhalf_pt.dQuat[j] * dt * 2;
                pt.Quat[j] = pt.Quat[j].normalize();

                pt.Vel[j] = vhalf_pt.Vel[j] + ptacc.acc[j] * dt;
                pt.AngSpd[j] = vhalf_pt.AngSpd[j] + ptacc.domg[j] * dt;

                double P_tmp = vhalf_pt.Vel[j].norm() * dt * 2;
                if (P_tmp > Vmax_thread[cur_th_num])Vmax_thread[cur_th_num] = P_tmp;//
            }
        }
        pt.Vmax = 0.0;//
        for (int j = 0; j < th_num; j++)
        {
            if (Vmax_thread[j] > pt.Vmax)
                pt.Vmax = Vmax_thread[j];
        }
        //
        for (int j = 0; j < bodyset.Num; j++)
        {
            bodyset.body[j].MassCenter += vhalf_bd.Vel[j] * dt * 2;
            bodyset.body[j].orien += vhalf_bd.dQuat[j] * dt * 2;
            bodyset.body[j].orien = bodyset.body[j].orien.normalize();

            bodyset.body[j].Vel = vhalf_bd.Vel[j] + bodyacc.acc[j] * dt;
            bodyset.body[j].AngularVel = vhalf_bd.AngSpd[j] + bodyacc.domg[j] * dt;

            double P_tmp_spin = bodyset.body[j].sphereRadius * (vhalf_bd.AngSpd[j].norm()) * dt * 2;
            double P_tmp_move = (vhalf_bd.Vel[j].norm()) * dt * 2;

            double P_tmp = P_tmp_spin + P_tmp_move;
            if (P_tmp > pt.Vmax) pt.Vmax = P_tmp;
        }

        pt.DRmax = pt.DRmax + pt.Vmax;//tree
        if (pt.DRmax > C)//update mesh
        {
            CreatTree(pt, bodyset, grid);
            pt.DRmax = 0.0;
        }

        bodyset.updateMatrix();
        force.ForceCaculate(
            time, pt, bodyset, control, mechPP, mechPT, mechPW, ode,
            part, poly_part, grid, vhalf_pt, vhalf_bd);

        for (int j = 0; j < pt.Num; j++)
        {
            DynEqnParticle(j, pt, force.FR[j], ptacc.acc[j], ptacc.domg[j]);
            vhalf_pt.Vel[j] += ptacc.acc[j] * 2*dt;
            vhalf_pt.AngSpd[j] += ptacc.domg[j] * 2*dt;
            vhalf_pt.quat[j] = (pt.Quat[j] + gama(pt.Quat[j], pt.AngSpd[j]) * dt).normalize();
            vhalf_pt.dQuat[j] = gama(vhalf_pt.quat[j], vhalf_pt.AngSpd[j]);
        }
        DynEqnBodySet(time, bodyset, force.poly_F, ddu);
        for (int j = 0; j < bodyset.Num; j++)
        {
            bodyacc.acc[j][0] = ddu[j][0]; bodyacc.acc[j][1] = ddu[j][1]; bodyacc.acc[j][2] = ddu[j][2];
            bodyacc.domg[j][0] = ddu[j][3]; bodyacc.domg[j][1] = ddu[j][4]; bodyacc.domg[j][2] = ddu[j][5];

            vhalf_bd.Vel[j] += bodyacc.acc[j] * dt * 2;
            vhalf_bd.AngSpd[j] += bodyacc.domg[j] * dt * 2;
            vhalf_bd.quat[j] = (bodyset.body[j].orien + gama(bodyset.body[j].orien, bodyset.body[j].AngularVel) * dt).normalize();
            vhalf_bd.dQuat[j] = gama(vhalf_bd.quat[j], vhalf_bd.AngSpd[j]);
        }

        // output check
        count++;
        if (count == ode.OutputInterval || i == ode.EndStep-1)//OutputInterval
        {
            cout << "  " << i + 1 << "/" << ode.EndStep << ",size = " << grid.size() << endl;

            if (control.OutputParticleForce) {
                force.ParticlesForceOutput(i, pt.Num, bodyset.Num);
            }
            if (control.OutputBodyForce) {
                force.BodysetForceOutput(i, pt.Num, bodyset.Num);
            }
            if (control.OutputParticleState) {
                pt.StateOutput(i);
            }
            if (control.OutputBodyState) {
                bodyset.StateOutput(time, force.poly_F, i);
            }
            
            count = 0;
        }
    }//ode.EndStep

    for (int i = 0; i < pt.Num; i++)
    {
        free(part[i]);
        free(poly_part[i]);
    }
    delete []part;
    delete []poly_part;
    for (int j = 0; j < bodyset.Num; j++) delete[]ddu[j];
    delete[]ddu;
    delete[] Vmax_thread;
}







#endif

void Simulation::summaryOutput(int time) {
	ofstream resfile;
	ostringstream convert;
	string Filename;

	Filename = "Data/DATA/OutputFile/Summary.par";

	const auto outputPath = file_utils::prepareOutputFile(
		Filename, "simulation summary output / 仿真摘要输出");
	resfile.open(outputPath); // 修改为默认模式（覆盖/清空重写），原为 ios::app (追加)
	if (!resfile.is_open()) {
		throw file_utils::pathError(
			"Cannot open / 无法打开",
			"simulation summary output / 仿真摘要输出",
			outputPath);
	}
	// 写入 UTF-8 BOM 以确保文本编辑器识别编码
	const char bom[] = { (char)0xEF, (char)0xBB, (char)0xBF };
	resfile.write(bom, sizeof(bom));
	
	//resfile << setiosflags(ios::scientific) << setprecision(PrecDouble); // output format
	//
	resfile <<"# Total time is "  << time <<" min"<< endl;
	resfile <<"# Particles' number is "  << pt.Num << endl;
	resfile <<"# Rigid bodies' number is " << bodyset.Num <<"\n"<< endl;
	resfile << "# Case directory is "
		<< file_utils::caseDirectory().generic_string() << endl;
	resfile <<endl;

	resfile << "# 2. Solving options setting (xxx format applied)" << endl;
	resfile << "ODE_Option         = "<< std::left<< setw(8) << ode.Option <<" # Solving option:  0- initial step checking-out; 1- do integration"<< endl;
	resfile << "ODE_StepSize       = "<< std::left<< setw(8) << ode.StepSize <<" # Step size of integration, unit: s"<< endl;
	resfile << "ODE_EndStep        = "<< std::left<< setw(8) << ode.EndStep <<" # Number of steps in intervals of StepSize"<< endl;
	resfile << "ODE_StartStep      = "<< std::left<< setw(8) << ode.StartStep <<" # Use this to change starting step numbering"<< endl;
	resfile << "ODE_OutputInterval = "<< std::left<< setw(8) << ode.OutputInterval <<" # Output interval\n"<< endl;
	resfile <<endl;

	resfile << "# 3. Mechanical parameters - 力学参数" << endl;
	resfile << "# PP: particle-particle / 颗粒-颗粒" << endl;
	writeMechanicalSet(resfile, "MECH_PP_", mechPP);
	resfile << "# PT: particle-triangle / 颗粒-三角面" << endl;
	writeMechanicalSet(resfile, "MECH_PT_", mechPT);
	resfile << "# PW: particle-wall / 颗粒-墙壁" << endl;
	writeMechanicalSet(resfile, "MECH_PW_", mechPW);
	resfile <<endl;

	resfile << "# Control Parameters - 控制参数" << endl;
	resfile << "CONTROL_GRAVITY_FLAG = " << std::left<< setw(8) << control.Grav_flag << " # 重力标志: 0-无, 1-常重力, 2-时变重力" << endl;
	resfile << "CONTROL_GRAVITY_X    = " << std::left<< setw(8) << control.g[0] << " # 重力X分量 / Gravity X component" << endl;
	resfile << "CONTROL_GRAVITY_Y    = " << std::left<< setw(8) << control.g[1] << " # 重力Y分量 / Gravity Y component" << endl;
	resfile << "CONTROL_GRAVITY_Z    = " << std::left<< setw(8) << control.g[2] << " # 重力Z分量 / Gravity Z component" << endl;
	resfile << "CONTROL_WALL_FLAG    = " << std::left<< setw(8) << control.Wall_flag << " # 墙壁标志" << endl;
	resfile << "CONTROL_MULTIBODY_FLAG = " << std::left<< setw(8) << control.Multibody_flag << " # 多体系统标志" << endl;
	resfile << "CONTROL_UNIVERSAL_GRAVITATION_FLAG = " << std::left<< setw(8) << control.Universal_gravitation_flag << " # 万有引力标志" << endl;
	resfile << "CONTROL_PARTICLES_FILE = " << control.ParticlesFile << " # Particles" << endl;
	resfile << "CONTROL_MULTIBODY_FILE = " << control.MultibodyFile << " # 多体系统文件 / Multibody system file" << endl;
	resfile << "CONTROL_WALL_FILE      = " << control.WallFile << " # 墙壁数据文件 / Wall data file" << endl;
	resfile << "CONTROL_OUTPUT_PARTICLE_STATE = " << control.OutputParticleState << endl;
	resfile << "CONTROL_OUTPUT_PARTICLE_FORCE = " << control.OutputParticleForce << endl;
	resfile << "CONTROL_OUTPUT_BODY_STATE     = " << control.OutputBodyState << endl;
	resfile << "CONTROL_OUTPUT_BODY_FORCE     = " << control.OutputBodyForce << endl;
	resfile << "CONTROL_OUTPUT_SUMMARY        = " << control.OutputSummary << endl;
	resfile <<endl;

	resfile.close();
}
