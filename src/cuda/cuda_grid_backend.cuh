#pragma once

#include "cuda_particle_backend.cuh"

#include <cuda_runtime.h>
#include <cstdint>

struct GpuGridArrays {
    int n = 0;
    double invMeshSize = 0.0;
    std::uint64_t* sortedCellKey = nullptr;
    int* sortedParticleId = nullptr;
    int* overflowCount = nullptr;
};

void gpuAllocateGrid(GpuGridArrays& grid, int particleCount);
void gpuFreeGrid(GpuGridArrays& grid);
void gpuBuildParticleGrid(GpuGridArrays& grid, const GpuParticleArrays& particles,
                          double meshSize);
