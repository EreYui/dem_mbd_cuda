#pragma once
#include <string>
#include "matrixtrans.h"
#include "particle.h"
#include "multibody.h"
#include "grid.h"
#include "wall.h"
#include "force.h"

struct CONTROL {
    int Grav_flag = 0;  //Grav_flag 0 zero gravity, 1 constant gravity, 2 gravity varying with time
    vector3d g = {0, 0, 0}; //gravity vector

    int Universal_gravitation_flag = 0; //Universal_gravitation_flag 0 no universal gravitation, 1 with universal gravitation

    //WallFile
    int Wall_flag = 0;
    std::string WallFile;
    WALL wall;

    //CouplePolyFile
    int Multibody_flag = 0;
    std::string MultibodyFile;
    std::string PrescribedMotionFile;

    //ParticlesFile
    std::string ParticlesFile;

    // Runtime output switches (0: disabled, 1: enabled).
    bool OutputParticleState = true;
    bool OutputParticleForce = true;
    bool OutputBodyState = true;
    bool OutputBodyForce = true;
    bool OutputSummary = true;
};


struct ODE {
    int Option;            // Solving option, 0: initial step checking-out; 1: do integration
    int EndStep;           //
    int StartStep;         // start step index
    int OutputInterval;    // output interval; file name like "res.000000000.bt" (9 digits)
    double StepSize;       //
};

struct MECH {
    double epsS = 1;               //tangential restitution coefficient
    double mu = 0.5;                 //tangential friction coefficient
    double epsN = 1;              //normal restitution coefficient
    double c = 1;                    //cohesive
    double beta = 1;               //
    double mu_T = 1.3;            //twist friction coefficient
    double mu_R = 1;           //rolling friction coefficient
    double kN = 1;
    double kS = 1;
    double cN = 1;
    double cS = 1;
    double kT = 1;
    double cT = 1;
    double kR = 1;
    double cR = 1;
};

/**
 * 仿真管理类
 */
class Simulation {
public:
    double time;

    PARTICLE pt;
    BODYSET bodyset;
    ODE ode;
    MECH mechPP; // particle-particle contact
    MECH mechPT; // particle-triangle contact
    MECH mechPW; // particle-wall contact
    CONTROL control;
    std::unordered_map<GridIndex, GridCell> grid;
    FORCE force;
    
    Simulation();
    ~Simulation();
    
    /**
     * 初始化仿真
     * @param configPath 配置文件路径 (可选)
     */
    void init(const std::string& parafile);
    
    /**
     * 运行仿真
     */
    void run();
    
    /**
     * 输出摘要信息
     */
    void summaryOutput(int time);

    void IntegrateDem();
    void IntegrateDemMultiBody();
private:
    //ODESolver solver;
    
    /**
     * 加载配置文件
     */
    void loadParas(const std::string& parafile);
};
