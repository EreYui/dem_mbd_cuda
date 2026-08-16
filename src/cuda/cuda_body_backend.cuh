#pragma once

#include <cuda_runtime.h>

struct GpuBodyStateArrays {
    int n = 0;
    double* cx = nullptr;
    double* cy = nullptr;
    double* cz = nullptr;
    double* qw = nullptr;
    double* qx = nullptr;
    double* qy = nullptr;
    double* qz = nullptr;
    double* vx = nullptr;
    double* vy = nullptr;
    double* vz = nullptr;
    double* wx = nullptr;
    double* wy = nullptr;
    double* wz = nullptr;
};

struct GpuTriangleArrays {
    int n = 0;
    int* bodyId = nullptr;
    int* triId = nullptr;

    double* lx1 = nullptr;
    double* ly1 = nullptr;
    double* lz1 = nullptr;
    double* lx2 = nullptr;
    double* ly2 = nullptr;
    double* lz2 = nullptr;
    double* lx3 = nullptr;
    double* ly3 = nullptr;
    double* lz3 = nullptr;
    double* lx4 = nullptr;
    double* ly4 = nullptr;
    double* lz4 = nullptr;
    double* lnx = nullptr;
    double* lny = nullptr;
    double* lnz = nullptr;

    double* wx1 = nullptr;
    double* wy1 = nullptr;
    double* wz1 = nullptr;
    double* wx2 = nullptr;
    double* wy2 = nullptr;
    double* wz2 = nullptr;
    double* wx3 = nullptr;
    double* wy3 = nullptr;
    double* wz3 = nullptr;
    double* wx4 = nullptr;
    double* wy4 = nullptr;
    double* wz4 = nullptr;
    double* wnx = nullptr;
    double* wny = nullptr;
    double* wnz = nullptr;
    double* aabbMinX = nullptr;
    double* aabbMinY = nullptr;
    double* aabbMinZ = nullptr;
    double* aabbMaxX = nullptr;
    double* aabbMaxY = nullptr;
    double* aabbMaxZ = nullptr;
};

void gpuAllocateBodyStates(GpuBodyStateArrays& b, int n);
void gpuFreeBodyStates(GpuBodyStateArrays& b);
void gpuAllocateTriangles(GpuTriangleArrays& t, int n);
void gpuFreeTriangles(GpuTriangleArrays& t);
void gpuTransformTriangles(GpuTriangleArrays& t, GpuBodyStateArrays& b);

