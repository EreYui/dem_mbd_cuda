#pragma once
#include "matrixtrans.h"


class WALL {
public:
    int number = 0;//
    vector3d* Orig = nullptr, * N = nullptr, * V = nullptr,
        * A = nullptr, * Omg = nullptr, * ON = nullptr;
//origin Orig, normal N, velocity of wall V(optional), oscillation amplitude A, oscillation angular frequency Omg, oscillation normal vector ON(x=Asin(Omg*t)ON)
    WALL() {}
    ~WALL() {}
    void LoadWalls(const std::string& parafile);

};
