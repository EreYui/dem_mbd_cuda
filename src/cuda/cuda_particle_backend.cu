#include "cuda_particle_backend.cuh"

#include <thrust/device_ptr.h>
#include <thrust/extrema.h>

#include <stdexcept>
#include <string>

namespace {

void checkCuda(cudaError_t err, const char* what)
{
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(err));
    }
}

void allocDouble(double*& ptr, int n)
{
    checkCuda(cudaMalloc(&ptr, sizeof(double) * n), "cudaMalloc(double)");
}

__global__ void clearForcesKernel(GpuParticleArrays d)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= d.n) return;
    d.fx[i] = 0.0;
    d.fy[i] = 0.0;
    d.fz[i] = 0.0;
    d.tx[i] = 0.0;
    d.ty[i] = 0.0;
    d.tz[i] = 0.0;
}

__global__ void addGravityKernel(GpuParticleArrays d, double gx, double gy, double gz)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= d.n) return;
    if (d.status[i] == -1) return;
    d.fx[i] += d.mass[i] * gx;
    d.fy[i] += d.mass[i] * gy;
    d.fz[i] += d.mass[i] * gz;
}

__global__ void addWallForcesKernel(GpuParticleArrays d, GpuWallArrays w, double kN, double epsN)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= d.n) return;
    if (d.status[i] == -1) return;

    const double lnepsN = log(epsN);
    const double denom = 9.86960440108935861883 + lnepsN * lnepsN;
    const double cN = -2.0 * lnepsN * sqrt(kN * d.mass[i] / denom);

    for (int iw = 0; iw < w.n; ++iw) {
        const double rx = d.x[i] - w.ox[iw];
        const double ry = d.y[i] - w.oy[iw];
        const double rz = d.z[i] - w.oz[iw];
        const double distance = rx * w.nx[iw] + ry * w.ny[iw] + rz * w.nz[iw];
        const double depth = d.radius[i] - distance;
        if (depth <= 0.0) continue;

        const double vn = d.vx[i] * w.nx[iw] + d.vy[i] * w.ny[iw] + d.vz[i] * w.nz[iw];
        const double fn = kN * depth - cN * vn;
        d.fx[i] += fn * w.nx[iw];
        d.fy[i] += fn * w.ny[iw];
        d.fz[i] += fn * w.nz[iw];
    }
}

__global__ void accelerationKernel(GpuParticleArrays d)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= d.n) return;
    d.ax[i] = d.fx[i] / d.mass[i];
    d.ay[i] = d.fy[i] / d.mass[i];
    d.az[i] = d.fz[i] / d.mass[i];
    const double qw = d.qw[i];
    const double qx = d.qx[i];
    const double qy = d.qy[i];
    const double qz = d.qz[i];
    // World-frame torque -> body-frame torque (Quat2DCM convention).
    const double tx = (1.0 - 2.0*(qy*qy + qz*qz))*d.tx[i]
                    + 2.0*(qx*qy + qw*qz)*d.ty[i]
                    + 2.0*(qx*qz - qw*qy)*d.tz[i];
    const double ty = 2.0*(qx*qy - qw*qz)*d.tx[i]
                    + (1.0 - 2.0*(qx*qx + qz*qz))*d.ty[i]
                    + 2.0*(qy*qz + qw*qx)*d.tz[i];
    const double tz = 2.0*(qx*qz + qw*qy)*d.tx[i]
                    + 2.0*(qy*qz - qw*qx)*d.ty[i]
                    + (1.0 - 2.0*(qx*qx + qy*qy))*d.tz[i];
    d.awx[i] = tx / d.inertia[i];
    d.awy[i] = ty / d.inertia[i];
    d.awz[i] = tz / d.inertia[i];
}

__device__ void quatDerivative(
    double qw,
    double qx,
    double qy,
    double qz,
    double wx,
    double wy,
    double wz,
    double& dqw,
    double& dqx,
    double& dqy,
    double& dqz)
{
    dqw = 0.5 * (-qx * wx - qy * wy - qz * wz);
    dqx = 0.5 * ( qw * wx + qy * wz - qz * wy);
    dqy = 0.5 * ( qw * wy - qx * wz + qz * wx);
    dqz = 0.5 * ( qw * wz + qx * wy - qy * wx);
}

__device__ void normalizeQuat(double& qw, double& qx, double& qy, double& qz)
{
    const double n2 = qw * qw + qx * qx + qy * qy + qz * qz;
    if (n2 <= 1.0e-24) {
        qw = 1.0;
        qx = 0.0;
        qy = 0.0;
        qz = 0.0;
        return;
    }
    const double inv = 1.0 / sqrt(n2);
    qw *= inv;
    qx *= inv;
    qy *= inv;
    qz *= inv;
}

__global__ void initializeLeapfrogKernel(GpuParticleArrays d, double halfDt)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= d.n) return;

    d.vhx[i] = d.vx[i] + d.ax[i] * halfDt;
    d.vhy[i] = d.vy[i] + d.ay[i] * halfDt;
    d.vhz[i] = d.vz[i] + d.az[i] * halfDt;
    d.whx[i] = d.wx[i] + d.awx[i] * halfDt;
    d.why[i] = d.wy[i] + d.awy[i] * halfDt;
    d.whz[i] = d.wz[i] + d.awz[i] * halfDt;

    double qhw = d.qw[i], qhx = d.qx[i], qhy = d.qy[i], qhz = d.qz[i];
    double dqw, dqx, dqy, dqz;
    quatDerivative(qhw, qhx, qhy, qhz, d.wx[i], d.wy[i], d.wz[i],
                   dqw, dqx, dqy, dqz);
    qhw += dqw * halfDt;
    qhx += dqx * halfDt;
    qhy += dqy * halfDt;
    qhz += dqz * halfDt;
    normalizeQuat(qhw, qhx, qhy, qhz);
    quatDerivative(
        qhw, qhx, qhy, qhz,
        d.whx[i], d.why[i], d.whz[i],
        d.hdqw[i], d.hdqx[i], d.hdqy[i], d.hdqz[i]);
}

__global__ void advancePositionKernel(GpuParticleArrays d, double fullDt)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= d.n) return;
    d.displacement[i] = 0.0;
    if (d.status[i] == -1) return;

    d.x[i] += d.vhx[i] * fullDt;
    d.y[i] += d.vhy[i] * fullDt;
    d.z[i] += d.vhz[i] * fullDt;
    d.displacement[i] = sqrt(d.vhx[i]*d.vhx[i] + d.vhy[i]*d.vhy[i] +
                             d.vhz[i]*d.vhz[i]) * fullDt;

    d.qw[i] += d.hdqw[i] * fullDt;
    d.qx[i] += d.hdqx[i] * fullDt;
    d.qy[i] += d.hdqy[i] * fullDt;
    d.qz[i] += d.hdqz[i] * fullDt;
    normalizeQuat(d.qw[i], d.qx[i], d.qy[i], d.qz[i]);

    d.vx[i] = d.vhx[i] + d.ax[i] * (0.5 * fullDt);
    d.vy[i] = d.vhy[i] + d.ay[i] * (0.5 * fullDt);
    d.vz[i] = d.vhz[i] + d.az[i] * (0.5 * fullDt);
    d.wx[i] = d.whx[i] + d.awx[i] * (0.5 * fullDt);
    d.wy[i] = d.why[i] + d.awy[i] * (0.5 * fullDt);
    d.wz[i] = d.whz[i] + d.awz[i] * (0.5 * fullDt);
}

__global__ void advanceHalfVelocityKernel(GpuParticleArrays d, double fullDt)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= d.n) return;

    d.vhx[i] += d.ax[i] * fullDt;
    d.vhy[i] += d.ay[i] * fullDt;
    d.vhz[i] += d.az[i] * fullDt;
    d.whx[i] += d.awx[i] * fullDt;
    d.why[i] += d.awy[i] * fullDt;
    d.whz[i] += d.awz[i] * fullDt;

    const double halfDt = 0.5 * fullDt;
    double qhw = d.qw[i], qhx = d.qx[i], qhy = d.qy[i], qhz = d.qz[i];
    double dqw, dqx, dqy, dqz;
    quatDerivative(qhw, qhx, qhy, qhz, d.wx[i], d.wy[i], d.wz[i],
                   dqw, dqx, dqy, dqz);
    qhw += dqw * halfDt;
    qhx += dqx * halfDt;
    qhy += dqy * halfDt;
    qhz += dqz * halfDt;
    normalizeQuat(qhw, qhx, qhy, qhz);
    quatDerivative(
        qhw, qhx, qhy, qhz,
        d.whx[i], d.why[i], d.whz[i],
        d.hdqw[i], d.hdqx[i], d.hdqy[i], d.hdqz[i]);
}

int gridSize(int n)
{
    return (n + 255) / 256;
}

} // namespace

void gpuAllocateParticles(GpuParticleArrays& d, int n)
{
    d.n = n;
    if (n <= 0) return;

    allocDouble(d.mass, n);
    allocDouble(d.inertia, n);
    allocDouble(d.radius, n);
    checkCuda(cudaMalloc(&d.id, sizeof(int) * n), "cudaMalloc(particle id)");
    checkCuda(cudaMalloc(&d.status, sizeof(int) * n), "cudaMalloc(status)");

    allocDouble(d.x, n);
    allocDouble(d.y, n);
    allocDouble(d.z, n);
    allocDouble(d.vx, n);
    allocDouble(d.vy, n);
    allocDouble(d.vz, n);
    allocDouble(d.wx, n);
    allocDouble(d.wy, n);
    allocDouble(d.wz, n);
    allocDouble(d.qw, n);
    allocDouble(d.qx, n);
    allocDouble(d.qy, n);
    allocDouble(d.qz, n);

    allocDouble(d.fx, n);
    allocDouble(d.fy, n);
    allocDouble(d.fz, n);
    allocDouble(d.tx, n);
    allocDouble(d.ty, n);
    allocDouble(d.tz, n);

    allocDouble(d.ax, n);
    allocDouble(d.ay, n);
    allocDouble(d.az, n);
    allocDouble(d.awx, n);
    allocDouble(d.awy, n);
    allocDouble(d.awz, n);

    allocDouble(d.vhx, n);
    allocDouble(d.vhy, n);
    allocDouble(d.vhz, n);
    allocDouble(d.whx, n);
    allocDouble(d.why, n);
    allocDouble(d.whz, n);
    allocDouble(d.hdqw, n);
    allocDouble(d.hdqx, n);
    allocDouble(d.hdqy, n);
    allocDouble(d.hdqz, n);
    allocDouble(d.displacement, n);
}

void gpuFreeParticles(GpuParticleArrays& d)
{
    cudaFree(d.mass);
    cudaFree(d.inertia);
    cudaFree(d.radius);
    cudaFree(d.id);
    cudaFree(d.status);
    cudaFree(d.x);
    cudaFree(d.y);
    cudaFree(d.z);
    cudaFree(d.vx);
    cudaFree(d.vy);
    cudaFree(d.vz);
    cudaFree(d.wx);
    cudaFree(d.wy);
    cudaFree(d.wz);
    cudaFree(d.qw);
    cudaFree(d.qx);
    cudaFree(d.qy);
    cudaFree(d.qz);
    cudaFree(d.fx);
    cudaFree(d.fy);
    cudaFree(d.fz);
    cudaFree(d.tx);
    cudaFree(d.ty);
    cudaFree(d.tz);
    cudaFree(d.ax);
    cudaFree(d.ay);
    cudaFree(d.az);
    cudaFree(d.awx);
    cudaFree(d.awy);
    cudaFree(d.awz);
    cudaFree(d.vhx);
    cudaFree(d.vhy);
    cudaFree(d.vhz);
    cudaFree(d.whx);
    cudaFree(d.why);
    cudaFree(d.whz);
    cudaFree(d.hdqw);
    cudaFree(d.hdqx);
    cudaFree(d.hdqy);
    cudaFree(d.hdqz);
    cudaFree(d.displacement);
    d = GpuParticleArrays{};
}

void gpuAllocateWalls(GpuWallArrays& w, int n)
{
    w.n = n;
    if (n <= 0) return;
    allocDouble(w.ox, n);
    allocDouble(w.oy, n);
    allocDouble(w.oz, n);
    allocDouble(w.nx, n);
    allocDouble(w.ny, n);
    allocDouble(w.nz, n);
}

void gpuFreeWalls(GpuWallArrays& w)
{
    cudaFree(w.ox);
    cudaFree(w.oy);
    cudaFree(w.oz);
    cudaFree(w.nx);
    cudaFree(w.ny);
    cudaFree(w.nz);
    w = GpuWallArrays{};
}

void gpuClearParticleForces(GpuParticleArrays& d)
{
    clearForcesKernel<<<gridSize(d.n), 256>>>(d);
    checkCuda(cudaGetLastError(), "clearForcesKernel");
}

void gpuAddGravity(GpuParticleArrays& d, double gx, double gy, double gz)
{
    addGravityKernel<<<gridSize(d.n), 256>>>(d, gx, gy, gz);
    checkCuda(cudaGetLastError(), "addGravityKernel");
}

void gpuAddWallForces(GpuParticleArrays& d, GpuWallArrays& w, double kN, double epsN)
{
    if (w.n <= 0) return;
    addWallForcesKernel<<<gridSize(d.n), 256>>>(d, w, kN, epsN);
    checkCuda(cudaGetLastError(), "addWallForcesKernel");
}

void gpuComputeParticleAcceleration(GpuParticleArrays& d)
{
    accelerationKernel<<<gridSize(d.n), 256>>>(d);
    checkCuda(cudaGetLastError(), "accelerationKernel");
}

void gpuInitializeParticleLeapfrog(GpuParticleArrays& d, double halfDt)
{
    initializeLeapfrogKernel<<<gridSize(d.n), 256>>>(d, halfDt);
    checkCuda(cudaGetLastError(), "initializeLeapfrogKernel");
}

void gpuAdvanceParticlePosition(GpuParticleArrays& d, double fullDt)
{
    advancePositionKernel<<<gridSize(d.n), 256>>>(d, fullDt);
    checkCuda(cudaGetLastError(), "advancePositionKernel");
}

void gpuAdvanceParticleHalfVelocity(GpuParticleArrays& d, double fullDt)
{
    advanceHalfVelocityKernel<<<gridSize(d.n), 256>>>(d, fullDt);
    checkCuda(cudaGetLastError(), "advanceHalfVelocityKernel");
}

double gpuMaxParticleDisplacement(GpuParticleArrays& d)
{
    if (d.n <= 0) return 0.0;
    thrust::device_ptr<double> begin(d.displacement);
    return *thrust::max_element(begin, begin + d.n);
}
