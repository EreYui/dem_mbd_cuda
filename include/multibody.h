#pragma once
#include <array>
#include <limits>
#include <string>
#include "matrixtrans.h"

// Axis-aligned bounding box used by the broad-phase grid.
struct AABB {
    vector3d min;
    vector3d max;

    AABB()
        : min(vector3d::Constant(std::numeric_limits<double>::infinity())),
          max(vector3d::Constant(-std::numeric_limits<double>::infinity())) {}

    void include(const vector3d& point) {
        min = min.min(point);
        max = max.max(point);
    }

    template <std::size_t N>
    void include(const std::array<vector3d, N>& points) {
        for (const auto& point : points) {
            include(point);
        }
    }
};

struct STLTetMesh {
    vector3d v1, v2, v3, v4, normal;
    double L;
} ;

enum BodyMotionState : int {
    BODY_DYNAMIC = 0,
    BODY_FIXED = 1,
    BODY_PRESCRIBED = 2
};

struct BODY {
    STLTetMesh* Shape = nullptr; // triangle mesh data
    AABB* eAABB = nullptr; // per-triangle AABB

    int pointNum = 0; // vertex count
    vector3d *point = nullptr; // vertices
    int Size = 0; // triangle count
    vector3i* tri = nullptr; // triangle indices
    double maxSize = 0;//
    double sphereRadius; // bounding sphere radius

    std::string name;
    // 0: force-driven; 1: fixed; 2: prescribed CSV or scheduled velocity.
    int state = BODY_DYNAMIC;
    double Mass = 1.0;
    vector3d MassCenter{ 0,0,0 };
    vector3d Vel{ 0,0,0 };
    vector3d AngularVel{ 0,0,0 }; // angular velocity
    vector4d orien{ 1,0,0,0 };
    matrix3d I = matrix3d::identity();
};


class BODYSET {
public:
    BODY* body = nullptr;
    int Num = 0;       // 刚体数量
    int SumFaceNum = 0;// 总面数
    matrix3d* IvDCM = nullptr;
    matrix3d* DCM = nullptr;

//
    vector3d* BodyLowerLmt = nullptr, * BodyUpperLmt = nullptr; // per-body AABB
    vector3d BodySetLowerLmt, BodySetUpperLmt; // overall AABB
    int*** triIdx = nullptr; // triangle adjacency: [body][tri][0..5]

    void updateMatrix() {
        for (int i = 0; i < Num; i++) {
            IvDCM[i] = Quat2IvDCM(body[i].orien);
            DCM[i] = Quat2DCM(body[i].orien);
        }
    }
    void StateOutput(double time, double** force, int i, double** impulse = nullptr);
    void LoadBodys(const std::string& filename);

private:
    
};
