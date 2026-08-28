#include "cuda/cuda_particle_backend.cuh"

#include <cuda_runtime.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void check(cudaError_t error, const char* message)
{
    if (error != cudaSuccess) {
        throw std::runtime_error(std::string(message) + ": " + cudaGetErrorString(error));
    }
}

template <typename T>
void upload(T* destination, const T& value)
{
    check(cudaMemcpy(destination, &value, sizeof(T), cudaMemcpyHostToDevice), "upload");
}

double download(const double* source)
{
    double value = 0.0;
    check(cudaMemcpy(&value, source, sizeof(double), cudaMemcpyDeviceToHost), "download");
    return value;
}

void requireClose(double observed, double expected, const char* name)
{
    if (std::abs(observed - expected) > 1.0e-14) {
        throw std::runtime_error(
            std::string(name) + " mismatch: observed=" + std::to_string(observed)
            + " expected=" + std::to_string(expected));
    }
}

} // namespace

int main()
{
    GpuParticleArrays owners;
    GpuComponentArrays components;
    try {
        gpuAllocateParticles(owners, 1);
        gpuAllocateComponents(components, 1, 1);
        const double halfSqrtTwo = std::sqrt(0.5);
        const double zero = 0.0;
        const double one = 1.0;
        const int owner = 0;
        upload(components.owner, owner);
        upload(components.bodyX, one);
        upload(components.bodyY, zero);
        upload(components.bodyZ, zero);
        upload(owners.qw, halfSqrtTwo);
        upload(owners.qx, zero);
        upload(owners.qy, zero);
        upload(owners.qz, halfSqrtTwo);
        const double px = 3.0, py = 4.0, pz = 5.0;
        upload(owners.x, px); upload(owners.y, py); upload(owners.z, pz);
        const double vx = 0.1, vy = 0.2, vz = 0.3;
        upload(owners.vx, vx); upload(owners.vy, vy); upload(owners.vz, vz);
        upload(owners.wx, zero); upload(owners.wy, one); upload(owners.wz, zero);
        const double vhx = 0.4, vhy = 0.5, vhz = 0.6;
        upload(owners.vhx, vhx); upload(owners.vhy, vhy); upload(owners.vhz, vhz);
        upload(owners.whx, one); upload(owners.why, zero); upload(owners.whz, zero);

        gpuUpdateComponents(components, owners);
        check(cudaDeviceSynchronize(), "component transform kernel");
        requireClose(download(components.offsetX), 0.0, "offset x");
        requireClose(download(components.offsetY), 1.0, "offset y");
        requireClose(download(components.offsetZ), 0.0, "offset z");
        requireClose(download(components.x), 3.0, "position x");
        requireClose(download(components.y), 5.0, "position y");
        requireClose(download(components.z), 5.0, "position z");
        // body omega (0,1,0) rotates to world (-1,0,0), so omega x r=(0,0,-1).
        requireClose(download(components.vx), 0.1, "velocity x");
        requireClose(download(components.vy), 0.2, "velocity y");
        requireClose(download(components.vz), -0.7, "velocity z");
        // half-step body omega (1,0,0) rotates parallel to the world offset.
        requireClose(download(components.vhx), 0.4, "half velocity x");
        requireClose(download(components.vhy), 0.5, "half velocity y");
        requireClose(download(components.vhz), 0.6, "half velocity z");
        gpuFreeComponents(components);
        gpuFreeParticles(owners);
    }
    catch (const std::exception& error) {
        gpuFreeComponents(components);
        gpuFreeParticles(owners);
        std::cerr << "component transform test failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    std::cout << "component transform tests passed\n";
    return EXIT_SUCCESS;
}
