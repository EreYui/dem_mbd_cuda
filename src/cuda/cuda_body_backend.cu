#include "cuda_body_backend.cuh"

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

void allocInt(int*& ptr, int n)
{
    checkCuda(cudaMalloc(&ptr, sizeof(int) * n), "cudaMalloc(int)");
}

__device__ void rotateByQuatIv(
    double qw,
    double qx,
    double qy,
    double qz,
    double x,
    double y,
    double z,
    double& ox,
    double& oy,
    double& oz)
{
    const double r00 = 1.0 - 2.0 * (qy * qy + qz * qz);
    const double r01 = 2.0 * (qx * qy - qw * qz);
    const double r02 = 2.0 * (qx * qz + qw * qy);
    const double r10 = 2.0 * (qx * qy + qw * qz);
    const double r11 = 1.0 - 2.0 * (qx * qx + qz * qz);
    const double r12 = 2.0 * (qy * qz - qw * qx);
    const double r20 = 2.0 * (qx * qz - qw * qy);
    const double r21 = 2.0 * (qy * qz + qw * qx);
    const double r22 = 1.0 - 2.0 * (qx * qx + qy * qy);

    ox = r00 * x + r01 * y + r02 * z;
    oy = r10 * x + r11 * y + r12 * z;
    oz = r20 * x + r21 * y + r22 * z;
}

__device__ void includePoint(
    double x,
    double y,
    double z,
    double& minX,
    double& minY,
    double& minZ,
    double& maxX,
    double& maxY,
    double& maxZ)
{
    minX = fmin(minX, x);
    minY = fmin(minY, y);
    minZ = fmin(minZ, z);
    maxX = fmax(maxX, x);
    maxY = fmax(maxY, y);
    maxZ = fmax(maxZ, z);
}

__global__ void transformTrianglesKernel(GpuTriangleArrays t, GpuBodyStateArrays b)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= t.n) return;

    const int ib = t.bodyId[i];
    const double qw = b.qw[ib];
    const double qx = b.qx[ib];
    const double qy = b.qy[ib];
    const double qz = b.qz[ib];

    double x1, y1, z1;
    double x2, y2, z2;
    double x3, y3, z3;
    double x4, y4, z4;
    double nx, ny, nz;

    rotateByQuatIv(qw, qx, qy, qz, t.lx1[i], t.ly1[i], t.lz1[i], x1, y1, z1);
    rotateByQuatIv(qw, qx, qy, qz, t.lx2[i], t.ly2[i], t.lz2[i], x2, y2, z2);
    rotateByQuatIv(qw, qx, qy, qz, t.lx3[i], t.ly3[i], t.lz3[i], x3, y3, z3);
    rotateByQuatIv(qw, qx, qy, qz, t.lx4[i], t.ly4[i], t.lz4[i], x4, y4, z4);
    rotateByQuatIv(qw, qx, qy, qz, t.lnx[i], t.lny[i], t.lnz[i], nx, ny, nz);

    x1 += b.cx[ib]; y1 += b.cy[ib]; z1 += b.cz[ib];
    x2 += b.cx[ib]; y2 += b.cy[ib]; z2 += b.cz[ib];
    x3 += b.cx[ib]; y3 += b.cy[ib]; z3 += b.cz[ib];
    x4 += b.cx[ib]; y4 += b.cy[ib]; z4 += b.cz[ib];

    t.wx1[i] = x1; t.wy1[i] = y1; t.wz1[i] = z1;
    t.wx2[i] = x2; t.wy2[i] = y2; t.wz2[i] = z2;
    t.wx3[i] = x3; t.wy3[i] = y3; t.wz3[i] = z3;
    t.wx4[i] = x4; t.wy4[i] = y4; t.wz4[i] = z4;
    t.wnx[i] = nx; t.wny[i] = ny; t.wnz[i] = nz;

    double minX = x1, minY = y1, minZ = z1;
    double maxX = x1, maxY = y1, maxZ = z1;
    includePoint(x2, y2, z2, minX, minY, minZ, maxX, maxY, maxZ);
    includePoint(x3, y3, z3, minX, minY, minZ, maxX, maxY, maxZ);
    includePoint(x4, y4, z4, minX, minY, minZ, maxX, maxY, maxZ);
    t.aabbMinX[i] = minX; t.aabbMinY[i] = minY; t.aabbMinZ[i] = minZ;
    t.aabbMaxX[i] = maxX; t.aabbMaxY[i] = maxY; t.aabbMaxZ[i] = maxZ;
}

int blockCount(int n)
{
    return (n + 255) / 256;
}

} // namespace

void gpuAllocateBodyStates(GpuBodyStateArrays& b, int n)
{
    b.n = n;
    if (n <= 0) return;
    allocDouble(b.cx, n); allocDouble(b.cy, n); allocDouble(b.cz, n);
    allocDouble(b.qw, n); allocDouble(b.qx, n); allocDouble(b.qy, n); allocDouble(b.qz, n);
    allocDouble(b.vx, n); allocDouble(b.vy, n); allocDouble(b.vz, n);
    allocDouble(b.wx, n); allocDouble(b.wy, n); allocDouble(b.wz, n);
}

void gpuFreeBodyStates(GpuBodyStateArrays& b)
{
    cudaFree(b.cx); cudaFree(b.cy); cudaFree(b.cz);
    cudaFree(b.qw); cudaFree(b.qx); cudaFree(b.qy); cudaFree(b.qz);
    cudaFree(b.vx); cudaFree(b.vy); cudaFree(b.vz);
    cudaFree(b.wx); cudaFree(b.wy); cudaFree(b.wz);
    b = GpuBodyStateArrays{};
}

void gpuAllocateTriangles(GpuTriangleArrays& t, int n)
{
    t.n = n;
    if (n <= 0) return;
    allocInt(t.bodyId, n); allocInt(t.triId, n);
    allocDouble(t.lx1, n); allocDouble(t.ly1, n); allocDouble(t.lz1, n);
    allocDouble(t.lx2, n); allocDouble(t.ly2, n); allocDouble(t.lz2, n);
    allocDouble(t.lx3, n); allocDouble(t.ly3, n); allocDouble(t.lz3, n);
    allocDouble(t.lx4, n); allocDouble(t.ly4, n); allocDouble(t.lz4, n);
    allocDouble(t.lnx, n); allocDouble(t.lny, n); allocDouble(t.lnz, n);
    allocDouble(t.wx1, n); allocDouble(t.wy1, n); allocDouble(t.wz1, n);
    allocDouble(t.wx2, n); allocDouble(t.wy2, n); allocDouble(t.wz2, n);
    allocDouble(t.wx3, n); allocDouble(t.wy3, n); allocDouble(t.wz3, n);
    allocDouble(t.wx4, n); allocDouble(t.wy4, n); allocDouble(t.wz4, n);
    allocDouble(t.wnx, n); allocDouble(t.wny, n); allocDouble(t.wnz, n);
    allocDouble(t.aabbMinX, n); allocDouble(t.aabbMinY, n); allocDouble(t.aabbMinZ, n);
    allocDouble(t.aabbMaxX, n); allocDouble(t.aabbMaxY, n); allocDouble(t.aabbMaxZ, n);
}

void gpuFreeTriangles(GpuTriangleArrays& t)
{
    cudaFree(t.bodyId); cudaFree(t.triId);
    cudaFree(t.lx1); cudaFree(t.ly1); cudaFree(t.lz1);
    cudaFree(t.lx2); cudaFree(t.ly2); cudaFree(t.lz2);
    cudaFree(t.lx3); cudaFree(t.ly3); cudaFree(t.lz3);
    cudaFree(t.lx4); cudaFree(t.ly4); cudaFree(t.lz4);
    cudaFree(t.lnx); cudaFree(t.lny); cudaFree(t.lnz);
    cudaFree(t.wx1); cudaFree(t.wy1); cudaFree(t.wz1);
    cudaFree(t.wx2); cudaFree(t.wy2); cudaFree(t.wz2);
    cudaFree(t.wx3); cudaFree(t.wy3); cudaFree(t.wz3);
    cudaFree(t.wx4); cudaFree(t.wy4); cudaFree(t.wz4);
    cudaFree(t.wnx); cudaFree(t.wny); cudaFree(t.wnz);
    cudaFree(t.aabbMinX); cudaFree(t.aabbMinY); cudaFree(t.aabbMinZ);
    cudaFree(t.aabbMaxX); cudaFree(t.aabbMaxY); cudaFree(t.aabbMaxZ);
    t = GpuTriangleArrays{};
}

void gpuTransformTriangles(GpuTriangleArrays& t, GpuBodyStateArrays& b)
{
    transformTrianglesKernel<<<blockCount(t.n), 256>>>(t, b);
    checkCuda(cudaGetLastError(), "transformTrianglesKernel");
}

