#pragma once

#include <cuda_runtime.h>

struct GpuParticleArrays {
    int n = 0;
    double* mass = nullptr;
    double* inertiaX = nullptr;
    double* inertiaY = nullptr;
    double* inertiaZ = nullptr;
    double* radius = nullptr;
    int* id = nullptr;
    int* status = nullptr;

    double* x = nullptr;
    double* y = nullptr;
    double* z = nullptr;
    double* vx = nullptr;
    double* vy = nullptr;
    double* vz = nullptr;
    double* wx = nullptr;
    double* wy = nullptr;
    double* wz = nullptr;
    double* qw = nullptr;
    double* qx = nullptr;
    double* qy = nullptr;
    double* qz = nullptr;

    double* fx = nullptr;
    double* fy = nullptr;
    double* fz = nullptr;
    double* tx = nullptr;
    double* ty = nullptr;
    double* tz = nullptr;

    double* ax = nullptr;
    double* ay = nullptr;
    double* az = nullptr;
    double* awx = nullptr;
    double* awy = nullptr;
    double* awz = nullptr;

    double* vhx = nullptr;
    double* vhy = nullptr;
    double* vhz = nullptr;
    double* whx = nullptr;
    double* why = nullptr;
    double* whz = nullptr;
    double* hdqw = nullptr;
    double* hdqx = nullptr;
    double* hdqy = nullptr;
    double* hdqz = nullptr;
    double* displacement = nullptr;
};

struct GpuComponentArrays {
    int n = 0;
    int ownerCount = 0;
    int* id = nullptr;
    int* owner = nullptr;
    int* ownerStart = nullptr; // [ownerCount + 1]
    double* radius = nullptr;
    double* bodyX = nullptr;
    double* bodyY = nullptr;
    double* bodyZ = nullptr;
    double* offsetX = nullptr;
    double* offsetY = nullptr;
    double* offsetZ = nullptr;
    double* x = nullptr;
    double* y = nullptr;
    double* z = nullptr;
    double* vx = nullptr;
    double* vy = nullptr;
    double* vz = nullptr;
    double* vhx = nullptr;
    double* vhy = nullptr;
    double* vhz = nullptr;
};

struct GpuWallArrays {
    int n = 0;
    double* ox = nullptr;
    double* oy = nullptr;
    double* oz = nullptr;
    double* nx = nullptr;
    double* ny = nullptr;
    double* nz = nullptr;
};

void gpuAllocateParticles(GpuParticleArrays& d, int n);
void gpuFreeParticles(GpuParticleArrays& d);
void gpuAllocateComponents(GpuComponentArrays& d, int componentCount, int ownerCount);
void gpuFreeComponents(GpuComponentArrays& d);
void gpuUpdateComponents(GpuComponentArrays& d, const GpuParticleArrays& owners);
void gpuAllocateWalls(GpuWallArrays& w, int n);
void gpuFreeWalls(GpuWallArrays& w);
void gpuClearParticleForces(GpuParticleArrays& d);
void gpuAddGravity(GpuParticleArrays& d, double gx, double gy, double gz);
void gpuAddWallForces(GpuParticleArrays& d, GpuWallArrays& w, double kN, double epsN);
void gpuComputeParticleAcceleration(GpuParticleArrays& d);
void gpuInitializeParticleLeapfrog(GpuParticleArrays& d, double halfDt);
void gpuAdvanceParticlePosition(GpuParticleArrays& d, double fullDt);
void gpuAdvanceParticleHalfVelocity(GpuParticleArrays& d, double fullDt);
double gpuMaxParticleDisplacement(GpuParticleArrays& d);
