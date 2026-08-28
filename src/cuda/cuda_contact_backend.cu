#include "cuda_contact_backend.cuh"

#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/sort.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <stdexcept>
#include <string>

namespace {

constexpr double kPi = 3.1415926535897932384626433832795;
constexpr double kGravityConstant = 6.67430e-11;

void checkCuda(cudaError_t err, const char* what)
{
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(err));
    }
}

void allocateDouble(double*& ptr, std::size_t count, const char* what)
{
    if (count == 0) return;
    checkCuda(cudaMalloc(&ptr, sizeof(double) * count), what);
}

void allocateInt(int*& ptr, std::size_t count, const char* what)
{
    if (count == 0) return;
    checkCuda(cudaMalloc(&ptr, sizeof(int) * count), what);
}

std::size_t checkedProduct(std::size_t lhs, std::size_t rhs,
                           std::size_t maximum, const char* what)
{
    if (lhs != 0 && rhs > maximum / lhs) {
        throw std::runtime_error(std::string(what) + " exceeds CUDA index capacity");
    }
    return lhs * rhs;
}

struct DVec3 {
    double x, y, z;
};

struct DVec2 {
    double x, y;
};

__device__ DVec3 makeVec(double x, double y, double z) { return {x, y, z}; }
__device__ DVec3 add(DVec3 a, DVec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
__device__ DVec3 sub(DVec3 a, DVec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
__device__ DVec3 mul(DVec3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
__device__ double dot(DVec3 a, DVec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
__device__ DVec3 cross(DVec3 a, DVec3 b)
{
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
__device__ double normSquared(DVec3 a) { return dot(a, a); }
__device__ double norm(DVec3 a) { return sqrt(normSquared(a)); }
__device__ DVec3 normalized(DVec3 a)
{
    const double length = norm(a);
    return length > 1.0e-14 ? mul(a, 1.0 / length) : makeVec(0.0, 0.0, 0.0);
}

__device__ DVec3 loadPosition(const GpuParticleArrays& p, int i)
{
    return {p.x[i], p.y[i], p.z[i]};
}

__device__ DVec3 loadComponentPosition(const GpuComponentArrays& c, int i)
{
    return {c.x[i], c.y[i], c.z[i]};
}

__device__ DVec3 loadComponentOffset(const GpuComponentArrays& c, int i)
{
    return {c.offsetX[i], c.offsetY[i], c.offsetZ[i]};
}

__device__ DVec3 loadComponentVelocity(const GpuComponentArrays& c, int i)
{
    return {c.vx[i], c.vy[i], c.vz[i]};
}

__device__ DVec3 loadComponentHalfStepVelocity(const GpuComponentArrays& c, int i)
{
    return {c.vhx[i], c.vhy[i], c.vhz[i]};
}

__device__ DVec3 loadAngularVelocity(const GpuParticleArrays& p, int i)
{
    return {p.wx[i], p.wy[i], p.wz[i]};
}

__device__ DVec3 loadHalfStepAngularVelocity(const GpuParticleArrays& p, int i)
{
    return {p.whx[i], p.why[i], p.whz[i]};
}

__device__ DVec3 rotateIv(double qw, double qx, double qy, double qz, DVec3 v)
{
    return {
        (1.0 - 2.0*(qy*qy + qz*qz))*v.x + 2.0*(qx*qy - qw*qz)*v.y + 2.0*(qx*qz + qw*qy)*v.z,
        2.0*(qx*qy + qw*qz)*v.x + (1.0 - 2.0*(qx*qx + qz*qz))*v.y + 2.0*(qy*qz - qw*qx)*v.z,
        2.0*(qx*qz - qw*qy)*v.x + 2.0*(qy*qz + qw*qx)*v.y + (1.0 - 2.0*(qx*qx + qy*qy))*v.z
    };
}

__device__ DVec3 particleAngularWorld(const GpuParticleArrays& p, int i)
{
    return rotateIv(p.qw[i], p.qx[i], p.qy[i], p.qz[i], loadAngularVelocity(p, i));
}

__device__ DVec3 particleHalfStepAngularWorld(const GpuParticleArrays& p, int i)
{
    return rotateIv(
        p.qw[i], p.qx[i], p.qy[i], p.qz[i], loadHalfStepAngularVelocity(p, i));
}

__device__ DVec3 bodyAngularWorld(const GpuBodyStateArrays& b, int i)
{
    return rotateIv(b.qw[i], b.qx[i], b.qy[i], b.qz[i],
                    makeVec(b.wx[i], b.wy[i], b.wz[i]));
}

__device__ __host__ std::uint64_t packCell(int x, int y, int z)
{
    constexpr int bias = 1 << 20;
    constexpr std::uint64_t mask = (1ull << 21) - 1ull;
    return ((static_cast<std::uint64_t>(x + bias) & mask) << 42) |
           ((static_cast<std::uint64_t>(y + bias) & mask) << 21) |
           (static_cast<std::uint64_t>(z + bias) & mask);
}

__device__ int lowerBound(const std::uint64_t* keys, int count, std::uint64_t key)
{
    int first = 0;
    int length = count;
    while (length > 0) {
        const int half = length >> 1;
        const int middle = first + half;
        if (keys[middle] < key) {
            first = middle + 1;
            length -= half + 1;
        } else {
            length = half;
        }
    }
    return first;
}

__device__ int upperBound(const std::uint64_t* keys, int count, std::uint64_t key)
{
    int first = 0;
    int length = count;
    while (length > 0) {
        const int half = length >> 1;
        const int middle = first + half;
        if (key < keys[middle]) {
            length = half;
        } else {
            first = middle + 1;
            length -= half + 1;
        }
    }
    return first;
}

__device__ double atomicMaxDouble(double* address, double value)
{
    auto* bits = reinterpret_cast<unsigned long long*>(address);
    unsigned long long old = *bits;
    while (__longlong_as_double(old) < value) {
        const unsigned long long assumed = old;
        old = atomicCAS(bits, assumed, __double_as_longlong(value));
        if (old == assumed) break;
    }
    return __longlong_as_double(old);
}

__device__ int acquireHistorySlot(std::uint64_t* keys, int* lastSeen,
                                  double* values, int particle, std::uint64_t key,
                                  int step, int slotCount, int* overflow,
                                  int* highWater)
{
    const int begin = particle * slotCount;
    int empty = -1;
    for (int slot = 0; slot < slotCount; ++slot) {
        const int index = begin + slot;
        if (keys[index] == key) {
            lastSeen[index] = step;
            atomicMax(highWater, slot + 1);
            return index;
        }
        if (keys[index] == 0 && empty < 0) empty = index;
    }
    if (empty < 0) {
        atomicAdd(overflow, 1);
        return -1;
    }
    keys[empty] = key;
    lastSeen[empty] = step;
    for (int k = 0; k < 9; ++k) values[empty * 9 + k] = 0.0;
    atomicMax(highWater, empty - begin + 1);
    return empty;
}

__device__ void expireHistory(std::uint64_t* keys, int* lastSeen, double* values,
                              int particle, int step, int slotCount)
{
    const int begin = particle * slotCount;
    for (int slot = 0; slot < slotCount; ++slot) {
        const int index = begin + slot;
        if (keys[index] != 0 && lastSeen[index] < step - 1) {
            keys[index] = 0;
            lastSeen[index] = -1;
            for (int k = 0; k < 9; ++k) values[index * 9 + k] = 0.0;
        }
    }
}

__device__ DVec3 historyLoad(double* values, int slot, int category)
{
    const int base = slot * 9 + category * 3;
    return {values[base], values[base + 1], values[base + 2]};
}

__device__ void historyStore(double* values, int slot, int category, DVec3 v)
{
    const int base = slot * 9 + category * 3;
    values[base] = v.x;
    values[base + 1] = v.y;
    values[base + 2] = v.z;
}

__device__ DVec3 transportTangentHistory(DVec3 previous, DVec3 normal)
{
    const double previousNorm = norm(previous);
    DVec3 projected = sub(previous, mul(normal, dot(previous, normal)));
    const double projectedNorm = norm(projected);
    if (previousNorm <= 1.0e-30 || projectedNorm <= 1.0e-30) {
        return makeVec(0.0, 0.0, 0.0);
    }
    // Rotate the stored elastic force/torque into the current contact plane
    // without introducing artificial relaxation when the normal changes.
    return mul(projected, previousNorm / projectedNorm);
}

__device__ DVec3 transportAxialHistory(DVec3 previous, DVec3 normal)
{
    const double previousNorm = norm(previous);
    if (previousNorm <= 1.0e-30) return makeVec(0.0, 0.0, 0.0);
    const double alignment = dot(previous, normal);
    if (fabs(alignment) <= 1.0e-30) return makeVec(0.0, 0.0, 0.0);
    return mul(normal, alignment < 0.0 ? -previousNorm : previousNorm);
}

__device__ DVec3 limitedHistoryResponseRates(
    DVec3 historyRate, DVec3 responseRate, double stiffness,
    double damping, double dt, DVec3 previous, double limit,
    double* values, int slot, int category, bool storeHistory = true)
{
    const DVec3 delta = mul(historyRate, dt);
    const DVec3 elasticIncrement = mul(delta, stiffness);
    const DVec3 elasticTrial = add(elasticIncrement, previous);
    const DVec3 trial = add(mul(responseRate, damping), elasticTrial);
    const double trialNorm = norm(elasticTrial);
    DVec3 response;
    if (trialNorm > limit) {
        const double deltaNorm = norm(delta);
        if (deltaNorm > 1.0e-14) {
            response = mul(elasticTrial, limit / trialNorm);
        } else {
            const double previousNorm = norm(previous);
            response = previousNorm > 1.0e-14
                ? mul(previous, limit / previousNorm)
                : makeVec(0.0, 0.0, 0.0);
        }
        if (storeHistory) historyStore(values, slot, category, response);
    } else {
        response = trial;
        if (storeHistory) {
            historyStore(values, slot, category, add(previous, elasticIncrement));
        }
    }
    return response;
}

__device__ DVec3 limitedHistoryResponse(DVec3 rate, double stiffness,
                                        double damping, double dt,
                                        DVec3 previous, double limit,
                                        double* values, int slot, int category,
                                        bool storeHistory = true)
{
    return limitedHistoryResponseRates(
        rate, rate, stiffness, damping, dt, previous, limit,
        values, slot, category, storeHistory);
}

__device__ double cross2(DVec2 a, DVec2 b) { return a.x*b.y - a.y*b.x; }
__device__ double dot2(DVec2 a, DVec2 b) { return a.x*b.x + a.y*b.y; }
__device__ DVec2 add2(DVec2 a, DVec2 b) { return {a.x+b.x, a.y+b.y}; }
__device__ DVec2 sub2(DVec2 a, DVec2 b) { return {a.x-b.x, a.y-b.y}; }
__device__ DVec2 mul2(DVec2 a, double s) { return {a.x*s, a.y*s}; }

__device__ double circleEdgeArea(DVec2 a, DVec2 b, double radius)
{
    const DVec2 d = sub2(b, a);
    const double aa = dot2(d, d);
    double cuts[4] = {0.0, 1.0, 0.0, 0.0};
    int cutCount = 2;
    if (aa > 1.0e-30) {
        const double bb = 2.0 * dot2(a, d);
        const double cc = dot2(a, a) - radius * radius;
        const double discriminant = bb * bb - 4.0 * aa * cc;
        if (discriminant > 0.0) {
            const double root = sqrt(discriminant);
            const double t0 = (-bb - root) / (2.0 * aa);
            const double t1 = (-bb + root) / (2.0 * aa);
            if (t0 > 0.0 && t0 < 1.0) cuts[cutCount++] = t0;
            if (t1 > 0.0 && t1 < 1.0) cuts[cutCount++] = t1;
        }
    }
    for (int i = 0; i < cutCount; ++i) {
        for (int j = i + 1; j < cutCount; ++j) {
            if (cuts[j] < cuts[i]) {
                const double tmp = cuts[i]; cuts[i] = cuts[j]; cuts[j] = tmp;
            }
        }
    }
    double area = 0.0;
    for (int i = 0; i + 1 < cutCount; ++i) {
        const double ta = cuts[i];
        const double tb = cuts[i + 1];
        const DVec2 u = add2(a, mul2(d, ta));
        const DVec2 v = add2(a, mul2(d, tb));
        const DVec2 midpoint = add2(a, mul2(d, 0.5 * (ta + tb)));
        if (dot2(midpoint, midpoint) <= radius * radius) {
            area += 0.5 * cross2(u, v);
        } else {
            area += 0.5 * radius * radius * atan2(cross2(u, v), dot2(u, v));
        }
    }
    return area;
}

__device__ double circleTriangleArea(DVec3 center, DVec3 normal,
                                     DVec3 a, DVec3 b, DVec3 c, double radius)
{
    DVec3 axis = fabs(normal.x) < 0.8 ? makeVec(1.0, 0.0, 0.0)
                                           : makeVec(0.0, 1.0, 0.0);
    DVec3 u = normalized(cross(normal, axis));
    DVec3 v = cross(normal, u);
    const DVec3 ra = sub(a, center);
    const DVec3 rb = sub(b, center);
    const DVec3 rc = sub(c, center);
    const DVec2 pa{dot(ra, u), dot(ra, v)};
    const DVec2 pb{dot(rb, u), dot(rb, v)};
    const DVec2 pc{dot(rc, u), dot(rc, v)};
    double area = circleEdgeArea(pa, pb, radius) +
                  circleEdgeArea(pb, pc, radius) +
                  circleEdgeArea(pc, pa, radius);
    area = fabs(area);
    return fmin(area, kPi * radius * radius);
}

__device__ double bodyTriangleAreaScale(const GpuComponentArrays& components,
                                        const GpuTriangleArrays& triangles,
                                        int component, int triangle)
{
    const DVec3 normal = makeVec(triangles.wnx[triangle], triangles.wny[triangle],
                                 triangles.wnz[triangle]);
    const DVec3 a = makeVec(triangles.wx2[triangle], triangles.wy2[triangle],
                            triangles.wz2[triangle]);
    const DVec3 b = makeVec(triangles.wx3[triangle], triangles.wy3[triangle],
                            triangles.wz3[triangle]);
    const DVec3 c = makeVec(triangles.wx4[triangle], triangles.wy4[triangle],
                            triangles.wz4[triangle]);
    const DVec3 position = loadComponentPosition(components, component);
    const double distance = dot(normal, sub(position, a));
    if (distance <= 0.0) return 0.0;
    const double circleRadiusSquared =
        components.radius[component] * components.radius[component]
        - distance * distance;
    if (circleRadiusSquared <= 1.0e-30) return 0.0;
    const double circleRadius = sqrt(circleRadiusSquared);
    const DVec3 circleCenter = sub(position, mul(normal, distance));
    const double overlapArea =
        circleTriangleArea(circleCenter, normal, a, b, c, circleRadius);
    if (overlapArea <= 0.0) return 0.0;
    const double areaScale =
        fmin(1.0, overlapArea / (kPi * circleRadiusSquared));
    return areaScale > 1.0e-12 ? areaScale : 0.0;
}

__global__ void countTriangleCellsKernel(GpuTriangleGrid grid,
                                         GpuTriangleArrays triangles)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= grid.triangleCount) return;
    const double minXd = floor((triangles.aabbMinX[i] - grid.padding) * grid.invMeshSize);
    const double minYd = floor((triangles.aabbMinY[i] - grid.padding) * grid.invMeshSize);
    const double minZd = floor((triangles.aabbMinZ[i] - grid.padding) * grid.invMeshSize);
    const double maxXd = floor((triangles.aabbMaxX[i] + grid.padding) * grid.invMeshSize);
    const double maxYd = floor((triangles.aabbMaxY[i] + grid.padding) * grid.invMeshSize);
    const double maxZd = floor((triangles.aabbMaxZ[i] + grid.padding) * grid.invMeshSize);
    constexpr int minCell = -(1 << 20);
    constexpr int maxCell = (1 << 20) - 1;
    if (!isfinite(minXd) || !isfinite(minYd) || !isfinite(minZd) ||
        !isfinite(maxXd) || !isfinite(maxYd) || !isfinite(maxZd) ||
        minXd < minCell || minYd < minCell || minZd < minCell ||
        maxXd > maxCell || maxYd > maxCell || maxZd > maxCell) {
        atomicAdd(grid.overflowCount, 1);
        grid.counts[i] = 0;
        return;
    }
    const int minX = static_cast<int>(minXd);
    const int minY = static_cast<int>(minYd);
    const int minZ = static_cast<int>(minZd);
    const int maxX = static_cast<int>(maxXd);
    const int maxY = static_cast<int>(maxYd);
    const int maxZ = static_cast<int>(maxZd);
    const long long nx = static_cast<long long>(maxX) - minX + 1;
    const long long ny = static_cast<long long>(maxY) - minY + 1;
    const long long nz = static_cast<long long>(maxZ) - minZ + 1;
    if (nx <= 0 || ny <= 0 || nz <= 0 || nx > LLONG_MAX / ny ||
        nx * ny > LLONG_MAX / nz) {
        atomicAdd(grid.overflowCount, 1);
        grid.counts[i] = 0;
        return;
    }
    const long long count64 = nx * ny * nz;
    if (count64 > INT_MAX) {
        atomicAdd(grid.overflowCount, 1);
        grid.counts[i] = 0;
        return;
    }
    grid.counts[i] = static_cast<int>(count64);
}

__global__ void fillTriangleCellsKernel(GpuTriangleGrid grid,
                                        GpuTriangleArrays triangles)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= grid.triangleCount) return;
    const int minX = static_cast<int>(floor((triangles.aabbMinX[i] - grid.padding) * grid.invMeshSize));
    const int minY = static_cast<int>(floor((triangles.aabbMinY[i] - grid.padding) * grid.invMeshSize));
    const int minZ = static_cast<int>(floor((triangles.aabbMinZ[i] - grid.padding) * grid.invMeshSize));
    const int maxX = static_cast<int>(floor((triangles.aabbMaxX[i] + grid.padding) * grid.invMeshSize));
    const int maxY = static_cast<int>(floor((triangles.aabbMaxY[i] + grid.padding) * grid.invMeshSize));
    const int maxZ = static_cast<int>(floor((triangles.aabbMaxZ[i] + grid.padding) * grid.invMeshSize));
    const int output = grid.offsets[i];
    int local = 0;
    for (int x = minX; x <= maxX && local < grid.counts[i]; ++x) {
        for (int y = minY; y <= maxY && local < grid.counts[i]; ++y) {
            for (int z = minZ; z <= maxZ && local < grid.counts[i]; ++z) {
                grid.sortedCellKey[output + local] = packCell(x, y, z);
                grid.sortedTriangleId[output + local] = i;
                ++local;
            }
        }
    }
}

__global__ void clearForcesKernel(GpuForceArrays f)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int particle6 = f.particleCount * 6;
    for (int index = i; index < particle6; index += blockDim.x * gridDim.x) {
        f.contact[index] = 0.0;
        f.field[index] = 0.0;
        f.wall[index] = 0.0;
        f.total[index] = 0.0;
    }
    if (f.bodyParticle) {
        const int bodyParticle6 = f.bodyCount * particle6;
        for (int index = i; index < bodyParticle6; index += blockDim.x * gridDim.x)
            f.bodyParticle[index] = 0.0;
    }
    if (f.bodyDetail) {
        const int bodyDetail19 = f.bodyCount * f.particleCount * 19;
        for (int index = i; index < bodyDetail19; index += blockDim.x * gridDim.x)
            f.bodyDetail[index] = 0.0;
    }
    const int body6 = f.bodyCount * 6;
    for (int index = i; index < body6; index += blockDim.x * gridDim.x)
        f.bodyResult[index] = 0.0;
    if (i == 0) {
        *f.particleContactCount = 0;
        *f.bodyContactCount = 0;
        *f.particleHistoryHighWater = 0;
        *f.bodyHistoryHighWater = 0;
        *f.bodyMaxDepth = 0.0;
        *f.bodyContactFz = 0.0;
    }
}

__device__ void addComponents(double* destination, int base, DVec3 force, DVec3 torque)
{
    destination[base] += force.x;
    destination[base + 1] += force.y;
    destination[base + 2] += force.z;
    destination[base + 3] += torque.x;
    destination[base + 4] += torque.y;
    destination[base + 5] += torque.z;
}

__device__ void particleContact(const GpuParticleArrays& p,
                                const GpuComponentArrays& components,
                                GpuContactHistory history,
                                GpuForceArrays forces,
                                const GpuMechanicalParams& mech,
                                int componentI, int componentJ, int step,
                                DVec3& forceSum, DVec3& torqueSum)
{
    const int i = components.owner[componentI];
    const int j = components.owner[componentJ];
    if (i == j) return;
    const DVec3 pi = loadComponentPosition(components, componentI);
    const DVec3 pj = loadComponentPosition(components, componentJ);
    const DVec3 delta = sub(pj, pi);
    const double distanceSquared = normSquared(delta);
    const double radiusI = components.radius[componentI];
    const double radiusJ = components.radius[componentJ];
    const double radiusSum = radiusI + radiusJ;
    if (distanceSquared >= radiusSum * radiusSum || distanceSquared <= 1.0e-28) return;

    const double distance = sqrt(distanceSquared);
    const double depth = radiusSum - distance;
    const DVec3 direction = mul(delta, 1.0 / distance);
    const DVec3 endpointVelocityI = loadComponentVelocity(components, componentI);
    const DVec3 endpointVelocityJ = loadComponentVelocity(components, componentJ);
    const DVec3 endpointRelativeVelocity =
        sub(endpointVelocityJ, endpointVelocityI);
    const DVec3 midpointVelocityI = mech.dt > 0.0 && p.status[i] != -1
        ? loadComponentHalfStepVelocity(components, componentI) : endpointVelocityI;
    const DVec3 midpointVelocityJ = mech.dt > 0.0 && p.status[j] != -1
        ? loadComponentHalfStepVelocity(components, componentJ) : endpointVelocityJ;
    const DVec3 midpointRelativeVelocity =
        sub(midpointVelocityJ, midpointVelocityI);
    const DVec3 endpointOmegaI = particleAngularWorld(p, i);
    const DVec3 endpointOmegaJ = particleAngularWorld(p, j);
    const DVec3 midpointOmegaI = mech.dt > 0.0 && p.status[i] != -1
        ? particleHalfStepAngularWorld(p, i) : endpointOmegaI;
    const DVec3 midpointOmegaJ = mech.dt > 0.0 && p.status[j] != -1
        ? particleHalfStepAngularWorld(p, j) : endpointOmegaJ;

    const double reducedMass = p.mass[i] * p.mass[j] / (p.mass[i] + p.mass[j]);
    const double kS = 2.0 * mech.kN / 7.0;
    const double logN = log(mech.epsN);
    const double cN = -2.0 * logN * sqrt(mech.kN * reducedMass /
                                         (kPi*kPi + logN*logN));
    const double logS = log(mech.epsS);
    const double cS = -2.0 * logS * sqrt(kS * reducedMass /
                                         (kPi*kPi + logS*logS));
    const double effectiveRadius = radiusI * radiusJ / radiusSum;
    const double betaRadiusSquared = mech.beta * effectiveRadius * mech.beta * effectiveRadius;
    const double kT = 2.0 * kS * betaRadiusSquared;
    const double cT = 2.0 * cS * betaRadiusSquared;
    const double kR = mech.kN * betaRadiusSquared;
    const double cR = cN * betaRadiusSquared;

    DVec3 normalForce = add(mul(direction, -mech.kN * depth),
                            mul(direction, cN * dot(endpointRelativeVelocity, direction)));
    normalForce = add(normalForce,
                      mul(direction, mech.cohesion * 4.0 * betaRadiusSquared));

    const std::uint64_t key = static_cast<std::uint64_t>(
        static_cast<std::uint32_t>(components.id[componentJ])) + 1ull;
    const int slot = acquireHistorySlot(history.particleKeys,
                                        history.particleLastSeen,
                                        history.particleValues,
                                        componentI, key, step,
                                        kGpuParticleHistorySlots,
                                        forces.particleHistoryOverflowCount,
                                        forces.particleHistoryHighWater);
    if (slot < 0) return;

    const double li = (radiusI*radiusI - radiusJ*radiusJ +
                       distanceSquared) / (2.0 * distance);
    const double lj = distance - li;
    DVec3 endpointTangentialVelocity = sub(
        endpointRelativeVelocity,
        mul(direction, dot(endpointRelativeVelocity, direction)));
    endpointTangentialVelocity = add(
        endpointTangentialVelocity, mul(cross(direction, endpointOmegaI), li));
    endpointTangentialVelocity = add(
        endpointTangentialVelocity, mul(cross(direction, endpointOmegaJ), lj));
    DVec3 midpointTangentialVelocity = sub(
        midpointRelativeVelocity,
        mul(direction, dot(midpointRelativeVelocity, direction)));
    midpointTangentialVelocity = add(
        midpointTangentialVelocity, mul(cross(direction, midpointOmegaI), li));
    midpointTangentialVelocity = add(
        midpointTangentialVelocity, mul(cross(direction, midpointOmegaJ), lj));

    const DVec3 oldTangent = transportTangentHistory(
        historyLoad(history.particleValues, slot, 1), direction);
    const DVec3 friction = limitedHistoryResponseRates(
        midpointTangentialVelocity, endpointTangentialVelocity,
        kS, cS, mech.dt, oldTangent,
        mech.mu * norm(normalForce), history.particleValues, slot, 1);

    const DVec3 endpointRelativeOmega = sub(endpointOmegaJ, endpointOmegaI);
    const DVec3 midpointRelativeOmega = sub(midpointOmegaJ, midpointOmegaI);
    const DVec3 endpointTwistRate =
        mul(direction, dot(endpointRelativeOmega, direction));
    const DVec3 midpointTwistRate =
        mul(direction, dot(midpointRelativeOmega, direction));
    const DVec3 endpointRollRate = sub(endpointRelativeOmega, endpointTwistRate);
    const DVec3 midpointRollRate = sub(midpointRelativeOmega, midpointTwistRate);
    const double momentScale = mech.beta * effectiveRadius * norm(normalForce);
    const DVec3 twist = limitedHistoryResponseRates(
        midpointTwistRate, endpointTwistRate, kT, cT, mech.dt,
        transportAxialHistory(
            historyLoad(history.particleValues, slot, 2), direction),
        mech.mu * mech.muT * momentScale,
        history.particleValues, slot, 2);
    const DVec3 roll = limitedHistoryResponseRates(
        midpointRollRate, endpointRollRate, kR, cR, mech.dt,
        transportTangentHistory(
            historyLoad(history.particleValues, slot, 0), direction),
        mech.muR * momentScale,
        history.particleValues, slot, 0);

    const DVec3 contactForce = add(normalForce, friction);
    const DVec3 componentCouple =
        add(add(twist, roll), mul(cross(direction, friction), li));
    forceSum = add(forceSum, contactForce);
    torqueSum = add(
        torqueSum,
        add(componentCouple,
            cross(loadComponentOffset(components, componentI), contactForce)));
    atomicAdd(forces.particleContactCount, 1);
}

__device__ void bodyContact(const GpuParticleArrays& p,
                            const GpuComponentArrays& components,
                            const GpuBodyStateArrays& bodies,
                            const GpuTriangleArrays& triangles,
                            GpuContactHistory history,
                            GpuForceArrays forces,
                            const GpuMechanicalParams& mech,
                            int component, int triangle, int historySlot,
                            DVec3& forceSum, DVec3& torqueSum)
{
    const int particle = components.owner[component];
    const int body = triangles.bodyId[triangle];
    DVec3 normal = makeVec(triangles.wnx[triangle], triangles.wny[triangle],
                           triangles.wnz[triangle]);
    const DVec3 a = makeVec(triangles.wx2[triangle], triangles.wy2[triangle],
                            triangles.wz2[triangle]);
    const DVec3 b = makeVec(triangles.wx3[triangle], triangles.wy3[triangle],
                            triangles.wz3[triangle]);
    const DVec3 c = makeVec(triangles.wx4[triangle], triangles.wy4[triangle],
                            triangles.wz4[triangle]);
    const DVec3 position = loadComponentPosition(components, component);
    const double distance = dot(normal, sub(position, a));
    if (distance <= 0.0) return;
    const double depth = components.radius[component] - distance;
    if (depth <= 0.0) return;
    const double circleRadiusSquared = components.radius[component]
                                       * components.radius[component] -
                                       distance*distance;
    if (circleRadiusSquared <= 1.0e-30) return;
    const double circleRadius = sqrt(circleRadiusSquared);
    const DVec3 circleCenter = sub(position, mul(normal, distance));
    const double overlapArea = circleTriangleArea(circleCenter, normal, a, b, c,
                                                  circleRadius);
    if (overlapArea <= 0.0) return;
    const double areaScale = fmin(1.0, overlapArea /
                                  (kPi * circleRadiusSquared));
    // Polygon clipping can leave tiny positive round-off slivers where the
    // legacy overlap routine reports no intersection. They contribute no
    // meaningful force but would consume contact-history slots indefinitely.
    if (areaScale <= 1.0e-12) return;

    normal = mul(normal, -1.0); // inward normal, matching the legacy model
    const DVec3 center = makeVec(bodies.cx[body], bodies.cy[body], bodies.cz[body]);
    const DVec3 arm = sub(circleCenter, center);
    const DVec3 omegaParticle = particleAngularWorld(p, particle);
    const DVec3 omegaBody = bodyAngularWorld(bodies, body);
    const DVec3 bodyVelocity = add(makeVec(bodies.vx[body], bodies.vy[body],
                                           bodies.vz[body]), cross(omegaBody, arm));
    const DVec3 relativeVelocity =
        sub(bodyVelocity, loadComponentVelocity(components, component));

    const double kS = 2.0 * mech.kN / 7.0;
    const double reducedMass = p.mass[particle];
    const double logN = log(mech.epsN);
    const double cN = -2.0 * logN * sqrt(mech.kN * reducedMass /
                                         (kPi*kPi + logN*logN));
    const double logS = log(mech.epsS);
    const double cS = -2.0 * logS * sqrt(kS * reducedMass /
                                         (kPi*kPi + logS*logS));
    const double radius = components.radius[component];
    const double betaRadiusSquared = mech.beta * radius * mech.beta * radius;
    const double kT = 2.0 * kS * betaRadiusSquared;
    const double cT = 2.0 * cS * betaRadiusSquared;
    const double kR = mech.kN * betaRadiusSquared;
    const double cR = cN * betaRadiusSquared;

    DVec3 normalForce = add(mul(normal, -mech.kN * depth),
                            mul(normal, cN * dot(relativeVelocity, normal)));
    normalForce = add(normalForce,
                      mul(normal, mech.cohesion * 4.0 * betaRadiusSquared));

    DVec3 tangentialVelocity = sub(relativeVelocity,
                                   mul(normal, dot(relativeVelocity, normal)));
    tangentialVelocity = add(tangentialVelocity,
                             mul(cross(normal, omegaParticle), distance));
    const DVec3 friction = limitedHistoryResponse(
        tangentialVelocity, kS, cS, mech.dt,
        transportTangentHistory(
            historyLoad(history.bodyValues, historySlot, 1), normal),
        mech.mu * norm(normalForce), history.bodyValues, historySlot, 1,
        false);

    const DVec3 relativeOmega = sub(omegaBody, omegaParticle);
    const DVec3 twistRate = mul(normal, dot(relativeOmega, normal));
    const DVec3 rollRate = add(sub(relativeOmega, twistRate),
                               mul(cross(normal, tangentialVelocity),
                                   0.5 / fmax(distance, 1.0e-14)));
    const double momentScale = mech.beta * radius * norm(normalForce);
    const DVec3 twist = limitedHistoryResponse(
        twistRate, kT, cT, mech.dt,
        transportAxialHistory(
            historyLoad(history.bodyValues, historySlot, 2), normal),
        mech.mu * mech.muT * momentScale,
        history.bodyValues, historySlot, 2, false);
    const DVec3 roll = limitedHistoryResponse(
        rollRate, kR, cR, mech.dt,
        transportTangentHistory(
            historyLoad(history.bodyValues, historySlot, 0), normal),
        mech.muR * momentScale,
        history.bodyValues, historySlot, 0, false);

    const DVec3 particleForce = mul(add(normalForce, friction), areaScale);
    const DVec3 contactMoment = add(twist, roll);
    const DVec3 componentCouple = mul(
        add(contactMoment, mul(cross(normal, friction), distance)), areaScale);
    const DVec3 particleTorque = add(
        componentCouple,
        cross(loadComponentOffset(components, component), particleForce));
    const DVec3 bodyNormalForce = mul(normalForce, -areaScale);
    const DVec3 bodyFriction = mul(friction, -areaScale);
    const DVec3 bodyNormalTorque = mul(cross(arm, normalForce), -areaScale);
    const DVec3 bodyFrictionTorque = mul(cross(arm, friction), -areaScale);
    const DVec3 bodyTwist = mul(twist, -areaScale);
    const DVec3 bodyRoll = mul(roll, -areaScale);
    const DVec3 bodyTorque = add(
        add(bodyNormalTorque, bodyFrictionTorque), add(bodyTwist, bodyRoll));

    forceSum = add(forceSum, particleForce);
    torqueSum = add(torqueSum, particleTorque);

    if (forces.bodyParticle) {
        const int bp = (body * forces.particleCount + particle) * 6;
        addComponents(forces.bodyParticle, bp, particleForce, particleTorque);
    }
    if (forces.bodyDetail) {
        const int detail = (body * forces.particleCount + particle) * 19;
        const DVec3 values[6] = {bodyNormalForce, bodyFriction, bodyNormalTorque,
                                 bodyFrictionTorque, bodyTwist, bodyRoll};
        for (int group = 0; group < 6; ++group) {
            forces.bodyDetail[detail + group*3] += values[group].x;
            forces.bodyDetail[detail + group*3 + 1] += values[group].y;
            forces.bodyDetail[detail + group*3 + 2] += values[group].z;
        }
        atomicMaxDouble(&forces.bodyDetail[detail + 18], depth);
    }

    atomicAdd(&forces.bodyResult[body*6], bodyNormalForce.x + bodyFriction.x);
    atomicAdd(&forces.bodyResult[body*6 + 1], bodyNormalForce.y + bodyFriction.y);
    atomicAdd(&forces.bodyResult[body*6 + 2], bodyNormalForce.z + bodyFriction.z);
    atomicAdd(&forces.bodyResult[body*6 + 3], bodyTorque.x);
    atomicAdd(&forces.bodyResult[body*6 + 4], bodyTorque.y);
    atomicAdd(&forces.bodyResult[body*6 + 5], bodyTorque.z);
    atomicAdd(forces.bodyContactCount, 1);
    atomicMaxDouble(forces.bodyMaxDepth, depth);
    atomicAdd(forces.bodyContactFz, bodyNormalForce.z + bodyFriction.z);
}

__device__ int advanceBodyManifoldHistory(
    const GpuParticleArrays& p, const GpuComponentArrays& components,
    const GpuBodyStateArrays& bodies,
    GpuContactHistory history, GpuForceArrays forces,
    const GpuMechanicalParams& mech, int component, int body, int step,
    DVec3 normal, DVec3 circleCenter, double distance)
{
    const int particle = components.owner[component];
    const std::uint64_t key = static_cast<std::uint64_t>(body + 1);
    const int slot = acquireHistorySlot(
        history.bodyKeys, history.bodyLastSeen, history.bodyValues,
        component, key, step, kGpuBodyHistorySlots,
        forces.bodyHistoryOverflowCount, forces.bodyHistoryHighWater);
    if (slot < 0) return -1;

    const DVec3 center =
        makeVec(bodies.cx[body], bodies.cy[body], bodies.cz[body]);
    const DVec3 arm = sub(circleCenter, center);
    const DVec3 endpointOmegaParticle = particleAngularWorld(p, particle);
    const bool dynamicPhysicalStep = mech.dt > 0.0 && p.status[particle] != -1;
    const DVec3 midpointOmegaParticle = dynamicPhysicalStep
        ? particleHalfStepAngularWorld(p, particle)
        : endpointOmegaParticle;
    const DVec3 omegaBody = bodyAngularWorld(bodies, body);
    const DVec3 bodyVelocity =
        add(makeVec(bodies.vx[body], bodies.vy[body], bodies.vz[body]),
            cross(omegaBody, arm));
    const DVec3 endpointParticleVelocity =
        loadComponentVelocity(components, component);
    const DVec3 midpointParticleVelocity = dynamicPhysicalStep
        ? loadComponentHalfStepVelocity(components, component)
        : endpointParticleVelocity;
    const DVec3 endpointRelativeVelocity =
        sub(bodyVelocity, endpointParticleVelocity);
    const DVec3 midpointRelativeVelocity =
        sub(bodyVelocity, midpointParticleVelocity);
    const double kS = 2.0 * mech.kN / 7.0;
    const double reducedMass = p.mass[particle];
    const double logN = log(mech.epsN);
    const double cN = -2.0 * logN * sqrt(
        mech.kN * reducedMass / (kPi*kPi + logN*logN));
    const double logS = log(mech.epsS);
    const double cS = -2.0 * logS * sqrt(
        kS * reducedMass / (kPi*kPi + logS*logS));
    const double radius = components.radius[component];
    const double betaRadiusSquared =
        mech.beta * radius * mech.beta * radius;
    const double kT = 2.0 * kS * betaRadiusSquared;
    const double cT = 2.0 * cS * betaRadiusSquared;
    const double kR = mech.kN * betaRadiusSquared;
    const double cR = cN * betaRadiusSquared;
    const double depth = radius - distance;
    DVec3 normalForce =
        add(mul(normal, -mech.kN * depth),
            mul(normal, cN * dot(endpointRelativeVelocity, normal)));
    normalForce = add(
        normalForce, mul(normal, mech.cohesion * 4.0 * betaRadiusSquared));

    DVec3 tangentialVelocity =
        sub(midpointRelativeVelocity,
            mul(normal, dot(midpointRelativeVelocity, normal)));
    tangentialVelocity = add(
        tangentialVelocity,
        mul(cross(normal, midpointOmegaParticle), distance));
    limitedHistoryResponse(
        tangentialVelocity, kS, cS, mech.dt,
        transportTangentHistory(
            historyLoad(history.bodyValues, slot, 1), normal),
        mech.mu * norm(normalForce), history.bodyValues, slot, 1, true);

    const DVec3 relativeOmega = sub(omegaBody, midpointOmegaParticle);
    const DVec3 twistRate =
        mul(normal, dot(relativeOmega, normal));
    const DVec3 rollRate =
        add(sub(relativeOmega, twistRate),
            mul(cross(normal, tangentialVelocity),
                0.5 / fmax(distance, 1.0e-14)));
    const double momentScale = mech.beta * radius * norm(normalForce);
    limitedHistoryResponse(
        twistRate, kT, cT, mech.dt,
        transportAxialHistory(
            historyLoad(history.bodyValues, slot, 2), normal),
        mech.mu * mech.muT * momentScale,
        history.bodyValues, slot, 2, true);
    limitedHistoryResponse(
        rollRate, kR, cR, mech.dt,
        transportTangentHistory(
            historyLoad(history.bodyValues, slot, 0), normal),
        mech.muR * momentScale, history.bodyValues, slot, 0, true);
    return slot;
}

__device__ DVec3 wallHistoryResponse(DVec3 rate, double stiffness,
                                     double damping, double dt,
                                     DVec3 previous, double limit,
                                     bool checkDampedTrial,
                                     bool zeroAtZeroIncrement,
                                     double* values, int slot, int category)
{
    const DVec3 delta = mul(rate, dt);
    const DVec3 elasticTrial = add(previous, mul(delta, stiffness));
    const DVec3 staticTrial = add(elasticTrial, mul(rate, damping));
    const double testNorm = norm(checkDampedTrial ? staticTrial : elasticTrial);
    if (testNorm > limit) {
        const double deltaNorm = norm(delta);
        const DVec3 limitedTrial = checkDampedTrial ? staticTrial : elasticTrial;
        const DVec3 response = deltaNorm > 1.0e-14 && testNorm > 1.0e-14
            ? mul(limitedTrial, limit / testNorm)
            : (zeroAtZeroIncrement ? makeVec(0.0, 0.0, 0.0) : previous);
        historyStore(values, slot, category, response);
        return response;
    }
    historyStore(values, slot, category, elasticTrial);
    return staticTrial;
}

__device__ void wallContact(const GpuParticleArrays& p,
                            const GpuComponentArrays& components,
                            const GpuWallArrays& walls,
                            GpuContactHistory history,
                            GpuForceArrays forces,
                            const GpuMechanicalParams& mech,
                            int component, int wall, int step,
                            DVec3& forceSum, DVec3& torqueSum)
{
    const int particle = components.owner[component];
    const DVec3 outward = makeVec(walls.nx[wall], walls.ny[wall], walls.nz[wall]);
    const DVec3 normal = mul(outward, -1.0);
    const DVec3 origin = makeVec(walls.ox[wall], walls.oy[wall], walls.oz[wall]);
    const DVec3 position = loadComponentPosition(components, component);
    const double distance = dot(sub(position, origin), outward);
    const double depth = components.radius[component] - distance;
    if (depth <= 0.0) return;

    const std::uint64_t key = static_cast<std::uint64_t>(wall) + 1ull;
    const int slot = acquireHistorySlot(
        history.wallKeys, history.wallLastSeen, history.wallValues,
        component, key, step, kGpuWallHistorySlots,
        forces.bodyHistoryOverflowCount, forces.bodyHistoryHighWater);
    if (slot < 0) return;

    const double radius = components.radius[component];
    const double lever = radius - depth;
    const DVec3 velocity = loadComponentVelocity(components, component);
    const DVec3 omega = particleAngularWorld(p, particle);
    const DVec3 relativeVelocity = mul(velocity, -1.0);
    DVec3 tangentialVelocity = sub(
        relativeVelocity, mul(normal, dot(relativeVelocity, normal)));
    tangentialVelocity = add(tangentialVelocity,
                             mul(cross(normal, omega), lever));
    const DVec3 relativeOmega = mul(omega, -1.0);
    const DVec3 twistRate = mul(normal, dot(relativeOmega, normal));
    DVec3 rollRate = sub(relativeOmega, twistRate);
    if (fabs(lever) > 1.0e-30) {
        rollRate = add(rollRate,
                       mul(cross(normal, tangentialVelocity), 0.5 / lever));
    }

    const double kS = 2.0 * mech.kN / 7.0;
    const double logN = log(mech.epsN);
    const double cN = -2.0 * logN * sqrt(
        mech.kN * p.mass[particle] / (kPi*kPi + logN*logN));
    const double logS = log(mech.epsS);
    const double cS = -2.0 * logS * sqrt(
        kS * p.mass[particle] / (kPi*kPi + logS*logS));
    const double betaRadiusSquared = mech.beta * radius * mech.beta * radius;
    const double kT = 2.0 * kS * betaRadiusSquared;
    const double cT = 2.0 * cS * betaRadiusSquared;
    const double kR = mech.kN * betaRadiusSquared;
    const double cR = cN * betaRadiusSquared;

    const DVec3 elasticNormal = mul(normal, -mech.kN * depth);
    const DVec3 dampingNormal = mul(
        normal, cN * dot(relativeVelocity, normal));
    const DVec3 normalForce = add(elasticNormal, dampingNormal);
    const double frictionLimit = mech.mu * norm(normalForce);
    const double momentScale = mech.beta * radius * norm(normalForce);

    const DVec3 friction = wallHistoryResponse(
        tangentialVelocity, kS, cS, mech.dt,
        transportTangentHistory(
            historyLoad(history.wallValues, slot, 1), normal), frictionLimit,
        false, false, history.wallValues, slot, 1);
    const DVec3 twist = wallHistoryResponse(
        twistRate, kT, cT, mech.dt,
        transportAxialHistory(
            historyLoad(history.wallValues, slot, 2), normal),
        mech.mu * mech.muT * momentScale,
        false, true, history.wallValues, slot, 2);
    const DVec3 roll = wallHistoryResponse(
        rollRate, kR, cR, mech.dt,
        transportTangentHistory(
            historyLoad(history.wallValues, slot, 0), normal),
        mech.muR * momentScale,
        true, true, history.wallValues, slot, 0);

    const DVec3 contactForce = add(normalForce, friction);
    const DVec3 componentCouple =
        add(add(twist, roll), mul(cross(normal, friction), lever));
    forceSum = add(forceSum, contactForce);
    torqueSum = add(
        torqueSum,
        add(componentCouple,
            cross(loadComponentOffset(components, component), contactForce)));
}

__global__ void computeContactsKernel(GpuParticleArrays p,
                                      GpuComponentArrays components,
                                      GpuGridArrays particleGrid,
                                      GpuWallArrays walls,
                                      GpuBodyStateArrays bodies,
                                      GpuTriangleArrays triangles,
                                      GpuTriangleGrid triangleGrid,
                                      GpuContactHistory history,
                                      GpuForceArrays forces,
                                      GpuMechanicalParams particleMech,
                                      GpuMechanicalParams bodyMech,
                                      GpuMechanicalParams wallMech,
                                      int step)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= p.n) return;

    DVec3 contactForce = makeVec(0.0, 0.0, 0.0);
    DVec3 contactTorque = makeVec(0.0, 0.0, 0.0);
    DVec3 fieldForce = makeVec(0.0, 0.0, 0.0);
    DVec3 wallForce = makeVec(0.0, 0.0, 0.0);
    DVec3 wallTorque = makeVec(0.0, 0.0, 0.0);
    DVec3 bodyForce = makeVec(0.0, 0.0, 0.0);
    DVec3 bodyTorque = makeVec(0.0, 0.0, 0.0);

    if (particleMech.gravityEnabled) {
        fieldForce = makeVec(p.mass[i]*particleMech.gx,
                             p.mass[i]*particleMech.gy,
                             p.mass[i]*particleMech.gz);
    }
    if (particleMech.universalGravityEnabled) {
        const DVec3 pi = loadPosition(p, i);
        for (int j = 0; j < p.n; ++j) {
            if (j == i) continue;
            const DVec3 delta = sub(pi, loadPosition(p, j));
            const double d2 = normSquared(delta);
            if (d2 > 1.0e-28) {
                fieldForce = add(fieldForce,
                    mul(delta, kGravityConstant*p.mass[i]*p.mass[j]/(d2*sqrt(d2))));
            }
        }
    }
    for (int component = components.ownerStart[i];
         component < components.ownerStart[i + 1]; ++component) {
        expireHistory(history.particleKeys, history.particleLastSeen,
                      history.particleValues, component, step,
                      kGpuParticleHistorySlots);
        expireHistory(history.wallKeys, history.wallLastSeen,
                      history.wallValues, component, step,
                      kGpuWallHistorySlots);
        expireHistory(history.bodyKeys, history.bodyLastSeen,
                      history.bodyValues, component, step,
                      kGpuBodyHistorySlots);

        for (int w = 0; w < walls.n; ++w) {
            wallContact(p, components, walls, history, forces, wallMech,
                        component, w, step, wallForce, wallTorque);
        }

        const DVec3 position = loadComponentPosition(components, component);
        const int ix = static_cast<int>(floor(position.x * particleGrid.invMeshSize));
        const int iy = static_cast<int>(floor(position.y * particleGrid.invMeshSize));
        const int iz = static_cast<int>(floor(position.z * particleGrid.invMeshSize));
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    const std::uint64_t key = packCell(ix+dx, iy+dy, iz+dz);
                    const int begin = lowerBound(particleGrid.sortedCellKey,
                                                 particleGrid.n, key);
                    const int end = upperBound(particleGrid.sortedCellKey,
                                               particleGrid.n, key);
                    for (int entry = begin; entry < end; ++entry) {
                        const int other = particleGrid.sortedEntityId[entry];
                        if (components.owner[other] != i) {
                            particleContact(
                                p, components, history, forces, particleMech,
                                component, other, step, contactForce, contactTorque);
                        }
                    }
                }
            }
        }

        if (triangles.n <= 0 || triangleGrid.entryCount <= 0) continue;
        const std::uint64_t key = packCell(ix, iy, iz);
        const int begin = lowerBound(
            triangleGrid.sortedCellKey, triangleGrid.entryCount, key);
        const int end = upperBound(
            triangleGrid.sortedCellKey, triangleGrid.entryCount, key);
        for (int body = 0; body < bodies.n; ++body) {
            double totalAreaScale = 0.0;
            DVec3 weightedNormal = makeVec(0.0, 0.0, 0.0);
            DVec3 weightedCircleCenter = makeVec(0.0, 0.0, 0.0);
            double weightedDistance = 0.0;
            for (int entry = begin; entry < end; ++entry) {
                const int triangle = triangleGrid.sortedTriangleId[entry];
                if (triangles.bodyId[triangle] != body) continue;
                const double areaScale = bodyTriangleAreaScale(
                    components, triangles, component, triangle);
                if (areaScale <= 0.0) continue;
                const DVec3 outward = makeVec(
                    triangles.wnx[triangle], triangles.wny[triangle],
                    triangles.wnz[triangle]);
                const DVec3 a = makeVec(
                    triangles.wx2[triangle], triangles.wy2[triangle],
                    triangles.wz2[triangle]);
                const double distance = dot(outward, sub(position, a));
                const DVec3 circleCenter = sub(position, mul(outward, distance));
                totalAreaScale += areaScale;
                weightedNormal = add(weightedNormal, mul(outward, -areaScale));
                weightedCircleCenter = add(
                    weightedCircleCenter, mul(circleCenter, areaScale));
                weightedDistance += distance * areaScale;
            }
            if (totalAreaScale <= 0.0) continue;
            const DVec3 manifoldNormal = normalized(weightedNormal);
            const DVec3 manifoldCircleCenter =
                mul(weightedCircleCenter, 1.0 / totalAreaScale);
            const double manifoldDistance = weightedDistance / totalAreaScale;
            const int historySlot = advanceBodyManifoldHistory(
                p, components, bodies, history, forces, bodyMech,
                component, body, step, manifoldNormal,
                manifoldCircleCenter, manifoldDistance);
            if (historySlot < 0) continue;

            GpuMechanicalParams readOnlyBodyMech = bodyMech;
            readOnlyBodyMech.dt = 0.0;
            for (int entry = begin; entry < end; ++entry) {
                const int triangle = triangleGrid.sortedTriangleId[entry];
                if (triangles.bodyId[triangle] != body) continue;
                bodyContact(
                    p, components, bodies, triangles, history, forces,
                    readOnlyBodyMech, component, triangle, historySlot,
                    bodyForce, bodyTorque);
            }
        }
    }

    const int base = i * 6;
    addComponents(forces.contact, base, contactForce, contactTorque);
    addComponents(forces.field, base, fieldForce, makeVec(0.0, 0.0, 0.0));
    addComponents(forces.wall, base, wallForce, wallTorque);
    const DVec3 totalForce = add(add(contactForce, fieldForce),
                                 add(wallForce, bodyForce));
    const DVec3 totalTorque = add(contactTorque, add(wallTorque, bodyTorque));
    addComponents(forces.total, base, totalForce, totalTorque);
    p.fx[i] = totalForce.x; p.fy[i] = totalForce.y; p.fz[i] = totalForce.z;
    p.tx[i] = totalTorque.x; p.ty[i] = totalTorque.y; p.tz[i] = totalTorque.z;
}

int blocks(int n) { return (n + 255) / 256; }

} // namespace

void gpuAllocateContactHistory(GpuContactHistory& h, int particleCount)
{
    if (particleCount <= 0 ||
        particleCount > INT_MAX / kGpuParticleHistorySlots ||
        particleCount > INT_MAX / kGpuWallHistorySlots ||
        particleCount > INT_MAX / kGpuBodyHistorySlots) {
        throw std::runtime_error("particle count exceeds CUDA contact-history index capacity");
    }
    h.particleCount = particleCount;
    const std::size_t particleSlots = static_cast<std::size_t>(particleCount) *
                                      kGpuParticleHistorySlots;
    const std::size_t wallSlots = static_cast<std::size_t>(particleCount) *
                                  kGpuWallHistorySlots;
    const std::size_t bodySlots = static_cast<std::size_t>(particleCount) *
                                  kGpuBodyHistorySlots;
    if (particleSlots == 0) return;
    checkCuda(cudaMalloc(&h.particleKeys, sizeof(std::uint64_t)*particleSlots),
              "cudaMalloc(particle history keys)");
    checkCuda(cudaMalloc(&h.wallKeys, sizeof(std::uint64_t)*wallSlots),
              "cudaMalloc(wall history keys)");
    checkCuda(cudaMalloc(&h.bodyKeys, sizeof(std::uint64_t)*bodySlots),
              "cudaMalloc(body history keys)");
    allocateInt(h.particleLastSeen, particleSlots, "cudaMalloc(particle history age)");
    allocateInt(h.wallLastSeen, wallSlots, "cudaMalloc(wall history age)");
    allocateInt(h.bodyLastSeen, bodySlots, "cudaMalloc(body history age)");
    allocateDouble(h.particleValues, particleSlots*9, "cudaMalloc(particle history values)");
    allocateDouble(h.wallValues, wallSlots*9, "cudaMalloc(wall history values)");
    allocateDouble(h.bodyValues, bodySlots*9, "cudaMalloc(body history values)");
    checkCuda(cudaMemset(h.particleKeys, 0, sizeof(std::uint64_t)*particleSlots),
              "clear particle history keys");
    checkCuda(cudaMemset(h.wallKeys, 0, sizeof(std::uint64_t)*wallSlots),
              "clear wall history keys");
    checkCuda(cudaMemset(h.bodyKeys, 0, sizeof(std::uint64_t)*bodySlots),
              "clear body history keys");
    checkCuda(cudaMemset(h.particleLastSeen, 0xff, sizeof(int)*particleSlots),
              "clear particle history age");
    checkCuda(cudaMemset(h.wallLastSeen, 0xff, sizeof(int)*wallSlots),
              "clear wall history age");
    checkCuda(cudaMemset(h.bodyLastSeen, 0xff, sizeof(int)*bodySlots),
              "clear body history age");
    checkCuda(cudaMemset(h.particleValues, 0, sizeof(double)*particleSlots*9),
              "clear particle history values");
    checkCuda(cudaMemset(h.wallValues, 0, sizeof(double)*wallSlots*9),
              "clear wall history values");
    checkCuda(cudaMemset(h.bodyValues, 0, sizeof(double)*bodySlots*9),
              "clear body history values");
}

void gpuFreeContactHistory(GpuContactHistory& h)
{
    cudaFree(h.particleKeys); cudaFree(h.wallKeys); cudaFree(h.bodyKeys);
    cudaFree(h.particleLastSeen); cudaFree(h.wallLastSeen); cudaFree(h.bodyLastSeen);
    cudaFree(h.particleValues); cudaFree(h.wallValues); cudaFree(h.bodyValues);
    h = GpuContactHistory{};
}

void gpuAllocateForces(GpuForceArrays& f, int particleCount, int bodyCount,
                       bool keepBodyParticle, bool keepBodyDetail)
{
    if (particleCount <= 0 || bodyCount < 0) {
        throw std::runtime_error("invalid particle/body count for CUDA force allocation");
    }
    f.particleCount = particleCount;
    f.bodyCount = bodyCount;
    const std::size_t kernelLimit = static_cast<std::size_t>(INT_MAX);
    const std::size_t particle6 = checkedProduct(
        static_cast<std::size_t>(particleCount), 6, kernelLimit,
        "particle force array");
    allocateDouble(f.contact, particle6, "cudaMalloc(contact forces)");
    allocateDouble(f.field, particle6, "cudaMalloc(field forces)");
    allocateDouble(f.wall, particle6, "cudaMalloc(wall forces)");
    allocateDouble(f.total, particle6, "cudaMalloc(total forces)");
    if (keepBodyParticle) {
        const std::size_t count = checkedProduct(
            static_cast<std::size_t>(bodyCount), particle6, kernelLimit,
            "body-particle output array");
        allocateDouble(f.bodyParticle,
                       count,
                       "cudaMalloc(body-particle forces)");
    }
    if (keepBodyDetail) {
        const std::size_t bodyParticles = checkedProduct(
            static_cast<std::size_t>(bodyCount),
            static_cast<std::size_t>(particleCount), kernelLimit,
            "body contact-detail output array");
        const std::size_t count = checkedProduct(
            bodyParticles, 19, kernelLimit,
            "body contact-detail output array");
        allocateDouble(f.bodyDetail,
                       count,
                       "cudaMalloc(body contact details)");
    }
    const std::size_t body6 = checkedProduct(
        static_cast<std::size_t>(bodyCount), 6, kernelLimit,
        "body resultant array");
    allocateDouble(f.bodyResult, body6,
                   "cudaMalloc(body resultant)");
    allocateInt(f.particleContactCount, 1, "cudaMalloc(particle contact count)");
    allocateInt(f.bodyContactCount, 1, "cudaMalloc(body contact count)");
    allocateInt(f.particleHistoryOverflowCount, 1,
                "cudaMalloc(particle history overflow count)");
    allocateInt(f.bodyHistoryOverflowCount, 1,
                "cudaMalloc(body history overflow count)");
    checkCuda(cudaMemset(f.particleHistoryOverflowCount, 0, sizeof(int)),
              "initialize particle history overflow count");
    checkCuda(cudaMemset(f.bodyHistoryOverflowCount, 0, sizeof(int)),
              "initialize body history overflow count");
    allocateInt(f.particleHistoryHighWater, 1,
                "cudaMalloc(particle history high-water mark)");
    allocateInt(f.bodyHistoryHighWater, 1,
                "cudaMalloc(body history high-water mark)");
    allocateDouble(f.bodyMaxDepth, 1, "cudaMalloc(body max depth)");
    allocateDouble(f.bodyContactFz, 1, "cudaMalloc(body contact Fz)");
}

void gpuFreeForces(GpuForceArrays& f)
{
    cudaFree(f.contact); cudaFree(f.field); cudaFree(f.wall); cudaFree(f.total);
    cudaFree(f.bodyParticle); cudaFree(f.bodyDetail); cudaFree(f.bodyResult);
    cudaFree(f.particleContactCount); cudaFree(f.bodyContactCount);
    cudaFree(f.particleHistoryOverflowCount);
    cudaFree(f.bodyHistoryOverflowCount);
    cudaFree(f.particleHistoryHighWater);
    cudaFree(f.bodyHistoryHighWater);
    cudaFree(f.bodyMaxDepth);
    cudaFree(f.bodyContactFz);
    f = GpuForceArrays{};
}

void gpuAllocateTriangleGrid(GpuTriangleGrid& g, int triangleCount)
{
    g.triangleCount = triangleCount;
    g.capacity = 0;
    if (triangleCount <= 0) return;
    allocateInt(g.counts, triangleCount, "cudaMalloc(triangle grid counts)");
    allocateInt(g.offsets, triangleCount, "cudaMalloc(triangle grid offsets)");
    allocateInt(g.overflowCount, 1, "cudaMalloc(triangle grid overflow)");
}

void gpuFreeTriangleGrid(GpuTriangleGrid& g)
{
    cudaFree(g.counts); cudaFree(g.offsets); cudaFree(g.sortedCellKey);
    cudaFree(g.sortedTriangleId); cudaFree(g.overflowCount);
    g = GpuTriangleGrid{};
}

void gpuBuildTriangleGrid(GpuTriangleGrid& g, const GpuTriangleArrays& triangles,
                          double meshSize, double padding)
{
    if (g.triangleCount <= 0) { g.entryCount = 0; return; }
    g.invMeshSize = 1.0 / meshSize;
    g.padding = padding;
    checkCuda(cudaMemset(g.overflowCount, 0, sizeof(int)),
              "clear triangle grid overflow");
    countTriangleCellsKernel<<<blocks(g.triangleCount), 256>>>(g, triangles);
    checkCuda(cudaGetLastError(), "countTriangleCellsKernel");
    thrust::device_ptr<int> counts(g.counts);
    thrust::device_ptr<int> offsets(g.offsets);
    const unsigned long long total = thrust::reduce(
        counts, counts + g.triangleCount, 0ull,
        thrust::plus<unsigned long long>());
    int overflow = 0;
    checkCuda(cudaMemcpy(&overflow, g.overflowCount, sizeof(int),
                         cudaMemcpyDeviceToHost),
              "copy triangle grid overflow");
    if (overflow != 0 || total > static_cast<unsigned long long>(INT_MAX)) {
        throw std::runtime_error(
            "triangle grid entry count exceeds the supported 32-bit index range");
    }
    g.entryCount = static_cast<int>(total);
    if (g.entryCount > g.capacity) {
        cudaFree(g.sortedCellKey);
        cudaFree(g.sortedTriangleId);
        g.sortedCellKey = nullptr;
        g.sortedTriangleId = nullptr;
        checkCuda(cudaMalloc(&g.sortedCellKey,
                             sizeof(std::uint64_t) *
                                 static_cast<std::size_t>(g.entryCount)),
                  "cudaMalloc(dynamic triangle grid keys)");
        allocateInt(g.sortedTriangleId, g.entryCount,
                    "cudaMalloc(dynamic triangle grid ids)");
        g.capacity = g.entryCount;
    }
    thrust::exclusive_scan(counts, counts + g.triangleCount, offsets);
    fillTriangleCellsKernel<<<blocks(g.triangleCount), 256>>>(g, triangles);
    checkCuda(cudaGetLastError(), "fillTriangleCellsKernel");
    thrust::device_ptr<std::uint64_t> keys(g.sortedCellKey);
    thrust::device_ptr<int> ids(g.sortedTriangleId);
    thrust::sort_by_key(keys, keys + g.entryCount, ids);
}

void gpuClearForces(GpuForceArrays& forces)
{
    clearForcesKernel<<<128, 256>>>(forces);
    checkCuda(cudaGetLastError(), "clearForcesKernel");
}

void gpuComputeContacts(const GpuParticleArrays& particles,
                        const GpuComponentArrays& components,
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
                        int stepIndex)
{
    if (particles.n <= 0) return;
    computeContactsKernel<<<blocks(particles.n), 256>>>(
        particles, components, particleGrid, walls, bodies, triangles, triangleGrid,
        history, forces, particleMech, bodyMech, wallMech, stepIndex);
    checkCuda(cudaGetLastError(), "computeContactsKernel");
}
