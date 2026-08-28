#include "cuda_grid_backend.cuh"

#include <thrust/device_ptr.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>

#include <stdexcept>
#include <string>

namespace {

void checkCuda(cudaError_t err, const char* what)
{
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(err));
    }
}

__device__ __host__ std::uint64_t packCell(int x, int y, int z)
{
    constexpr int bias = 1 << 20;
    constexpr std::uint64_t mask = (1ull << 21) - 1ull;
    return ((static_cast<std::uint64_t>(x + bias) & mask) << 42) |
           ((static_cast<std::uint64_t>(y + bias) & mask) << 21) |
           (static_cast<std::uint64_t>(z + bias) & mask);
}

__global__ void computeKeysKernel(GpuGridArrays grid, GpuParticleArrays p)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= p.n) return;
    constexpr int minCell = -(1 << 20);
    constexpr int maxCell = (1 << 20) - 1;
    const double cellX = floor(p.x[i] * grid.invMeshSize);
    const double cellY = floor(p.y[i] * grid.invMeshSize);
    const double cellZ = floor(p.z[i] * grid.invMeshSize);
    if (!isfinite(cellX) || !isfinite(cellY) || !isfinite(cellZ) ||
        cellX < minCell || cellX > maxCell ||
        cellY < minCell || cellY > maxCell ||
        cellZ < minCell || cellZ > maxCell) {
        atomicAdd(grid.overflowCount, 1);
        grid.sortedCellKey[i] = 0;
        return;
    }
    const int ix = static_cast<int>(cellX);
    const int iy = static_cast<int>(cellY);
    const int iz = static_cast<int>(cellZ);
    grid.sortedCellKey[i] = packCell(ix, iy, iz);
}

__global__ void computeComponentKeysKernel(
    GpuGridArrays grid, GpuComponentArrays components)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= components.n) return;
    constexpr int minCell = -(1 << 20);
    constexpr int maxCell = (1 << 20) - 1;
    const double cellX = floor(components.x[i] * grid.invMeshSize);
    const double cellY = floor(components.y[i] * grid.invMeshSize);
    const double cellZ = floor(components.z[i] * grid.invMeshSize);
    if (!isfinite(cellX) || !isfinite(cellY) || !isfinite(cellZ)
        || cellX < minCell || cellX > maxCell
        || cellY < minCell || cellY > maxCell
        || cellZ < minCell || cellZ > maxCell) {
        atomicAdd(grid.overflowCount, 1);
        grid.sortedCellKey[i] = 0;
        return;
    }
    grid.sortedCellKey[i] = packCell(
        static_cast<int>(cellX), static_cast<int>(cellY), static_cast<int>(cellZ));
}

int blocks(int n) { return (n + 255) / 256; }

} // namespace

void gpuAllocateGrid(GpuGridArrays& grid, int particleCount)
{
    grid.n = particleCount;
    if (particleCount <= 0) return;
    checkCuda(cudaMalloc(&grid.sortedCellKey,
                         sizeof(std::uint64_t) * particleCount),
              "cudaMalloc(particle grid keys)");
    checkCuda(cudaMalloc(&grid.sortedEntityId, sizeof(int) * particleCount),
              "cudaMalloc(particle grid ids)");
    checkCuda(cudaMalloc(&grid.overflowCount, sizeof(int)),
              "cudaMalloc(particle grid overflow counter)");
}

void gpuFreeGrid(GpuGridArrays& grid)
{
    cudaFree(grid.sortedCellKey);
    cudaFree(grid.sortedEntityId);
    cudaFree(grid.overflowCount);
    grid = GpuGridArrays{};
}

void gpuBuildParticleGrid(GpuGridArrays& grid, const GpuParticleArrays& particles,
                          double meshSize)
{
    if (grid.n <= 0) return;
    grid.invMeshSize = 1.0 / meshSize;
    thrust::device_ptr<std::uint64_t> keys(grid.sortedCellKey);
    thrust::device_ptr<int> ids(grid.sortedEntityId);
    thrust::sequence(ids, ids + grid.n);
    checkCuda(cudaMemset(grid.overflowCount, 0, sizeof(int)),
              "clear particle grid overflow counter");
    computeKeysKernel<<<blocks(grid.n), 256>>>(grid, particles);
    checkCuda(cudaGetLastError(), "computeKeysKernel");
    int overflow = 0;
    checkCuda(cudaMemcpy(&overflow, grid.overflowCount, sizeof(int),
                         cudaMemcpyDeviceToHost),
              "copy particle grid overflow counter");
    if (overflow != 0) {
        throw std::runtime_error(
            "particle grid coordinate exceeds signed 21-bit hash range");
    }
    thrust::sort_by_key(keys, keys + grid.n, ids);
}

void gpuBuildComponentGrid(GpuGridArrays& grid,
                           const GpuComponentArrays& components,
                           double meshSize)
{
    if (grid.n != components.n || grid.n <= 0 || meshSize <= 0.0) {
        throw std::runtime_error("invalid component grid shape or mesh size");
    }
    grid.invMeshSize = 1.0 / meshSize;
    thrust::device_ptr<std::uint64_t> keys(grid.sortedCellKey);
    thrust::device_ptr<int> ids(grid.sortedEntityId);
    thrust::sequence(ids, ids + grid.n);
    checkCuda(cudaMemset(grid.overflowCount, 0, sizeof(int)),
              "clear component grid overflow counter");
    computeComponentKeysKernel<<<blocks(grid.n), 256>>>(grid, components);
    checkCuda(cudaGetLastError(), "computeComponentKeysKernel");
    int overflow = 0;
    checkCuda(cudaMemcpy(&overflow, grid.overflowCount, sizeof(int),
                         cudaMemcpyDeviceToHost),
              "copy component grid overflow counter");
    if (overflow != 0) {
        throw std::runtime_error(
            "component grid coordinate exceeds signed 21-bit hash range");
    }
    thrust::sort_by_key(keys, keys + grid.n, ids);
}
