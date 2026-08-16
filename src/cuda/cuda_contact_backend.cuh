#pragma once

#include "cuda_body_backend.cuh"
#include "cuda_grid_backend.cuh"
#include "cuda_particle_backend.cuh"

#include <cuda_runtime.h>
#include <cstddef>
#include <cstdint>

constexpr int kGpuParticleHistorySlots = 64;
constexpr int kGpuWallHistorySlots = 16;
constexpr int kGpuBodyHistorySlots = 64;
constexpr int kGpuMaxCellsPerTriangle = 125;

struct GpuMechanicalParams {
    double epsS = 1.0;
    double mu = 0.5;
    double epsN = 1.0;
    double cohesion = 0.0;
    double beta = 1.0;
    double muT = 1.3;
    double muR = 1.0;
    double kN = 1.0;
    double dt = 0.0;
    double gx = 0.0;
    double gy = 0.0;
    double gz = 0.0;
    int gravityEnabled = 0;
    int universalGravityEnabled = 0;
};

struct GpuContactHistory {
    int particleCount = 0;
    std::uint64_t* particleKeys = nullptr;
    std::uint64_t* wallKeys = nullptr;
    std::uint64_t* bodyKeys = nullptr;
    int* particleLastSeen = nullptr;
    int* wallLastSeen = nullptr;
    int* bodyLastSeen = nullptr;
    double* particleValues = nullptr; // [particle][slot][roll,tangent,twist][xyz]
    double* wallValues = nullptr;
    double* bodyValues = nullptr;
};

struct GpuForceArrays {
    int particleCount = 0;
    int bodyCount = 0;
    double* contact = nullptr;      // [particle][6]
    double* field = nullptr;        // [particle][6]
    double* wall = nullptr;         // [particle][6]
    double* total = nullptr;        // [particle][6]
    double* bodyParticle = nullptr; // [body][particle][6], body -> particle
    double* bodyDetail = nullptr;   // [body][particle][19], force on body
    double* bodyResult = nullptr;   // [body][6], force on body
    int* particleContactCount = nullptr;
    int* bodyContactCount = nullptr;
    int* particleHistoryOverflowCount = nullptr;
    int* bodyHistoryOverflowCount = nullptr;
    int* particleHistoryHighWater = nullptr;
    int* bodyHistoryHighWater = nullptr;
    double* bodyMaxDepth = nullptr;
    double* bodyContactFz = nullptr;
};

struct GpuTriangleGrid {
    int triangleCount = 0;
    int capacity = 0;
    int entryCount = 0;
    double invMeshSize = 0.0;
    double padding = 0.0;
    int* counts = nullptr;
    int* offsets = nullptr;
    std::uint64_t* sortedCellKey = nullptr;
    int* sortedTriangleId = nullptr;
    int* overflowCount = nullptr;
};

void gpuAllocateContactHistory(GpuContactHistory& history, int particleCount);
void gpuFreeContactHistory(GpuContactHistory& history);
void gpuAllocateForces(GpuForceArrays& forces, int particleCount, int bodyCount,
                       bool keepBodyParticle, bool keepBodyDetail);
void gpuFreeForces(GpuForceArrays& forces);
void gpuAllocateTriangleGrid(GpuTriangleGrid& grid, int triangleCount);
void gpuFreeTriangleGrid(GpuTriangleGrid& grid);
void gpuBuildTriangleGrid(GpuTriangleGrid& grid, const GpuTriangleArrays& triangles,
                          double meshSize, double padding);
void gpuClearForces(GpuForceArrays& forces);
void gpuComputeContacts(const GpuParticleArrays& particles,
                        const GpuGridArrays& particleGrid,
                        const GpuWallArrays& walls,
                        const GpuBodyStateArrays& bodies,
                        const GpuTriangleArrays& triangles,
                        const GpuTriangleGrid& triangleGrid,
                        GpuContactHistory& history,
                        GpuForceArrays& forces,
                        const GpuMechanicalParams& particleMech,
                        const GpuMechanicalParams& bodyMech,
                        const GpuMechanicalParams& wallMech,
                        int stepIndex);
