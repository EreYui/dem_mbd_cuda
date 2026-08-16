#pragma once
#include <string>
#include "matrixtrans.h"

class PARTICLE {
public:
    int Num = 0;   // 颗粒数量
    double* Mass = nullptr, * Radius = nullptr, * Inertia = nullptr;// 质量, 半径, 转动惯量

    int* Number = nullptr;  //颗粒编号
    int *Status = nullptr;  //0:正常, 1:接触, 2:重叠
    vector3d* Pos = nullptr, * Vel = nullptr, * AngSpd = nullptr;//位置, 速度, 角速度
    vector4d* Quat = nullptr; // 四元数表示的姿态

    double maxRadius;//最大半径
    double Vmax;//?
    double MEmax; // max energy metric (legacy)
    double DRmax; // max displacement metric (legacy)

//
    double MeshSize; // mesh size
    vector3i MeshNum; // mesh dimensions

    vector3i* Idx = nullptr; // cell indices per particle
    vector3d LowerLmt, UpperLmt;

    void LoadParticles(std::string iniconfile);
    void StateOutput(int i);
    void FreeParticles();
    ~PARTICLE();
} ;
