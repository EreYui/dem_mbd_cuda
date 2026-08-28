#include "cuda_backend.h"

#include "cuda/cuda_body_backend.cuh"
#include "cuda/cuda_contact_backend.cuh"
#include "cuda/cuda_grid_backend.cuh"
#include "cuda/cuda_particle_backend.cuh"
#include "force.h"
#include "particle.h"
#include "simulation.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

void checkCuda(cudaError_t err, const char* what)
{
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(err));
    }
}

template <typename T>
void upload(T* destination, const T* source, std::size_t count, const char* what)
{
    if (count == 0) return;
    checkCuda(cudaMemcpy(destination, source, sizeof(T)*count,
                         cudaMemcpyHostToDevice), what);
}

template <typename T>
void download(T* destination, const T* source, std::size_t count,
              const char* what)
{
    if (count == 0) return;
    checkCuda(cudaMemcpy(destination, source, sizeof(T)*count,
                         cudaMemcpyDeviceToHost), what);
}

int countTriangles(const BODYSET& bodies)
{
    std::size_t total = 0;
    for (int body = 0; body < bodies.Num; ++body) {
        if (bodies.body[body].Size < 0) {
            throw std::runtime_error("body triangle count cannot be negative");
        }
        total += static_cast<std::size_t>(bodies.body[body].Size);
        if (total > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("total triangle count exceeds CUDA 32-bit index capacity");
        }
    }
    return static_cast<int>(total);
}

void validateStableParticleIds(const PARTICLE& p)
{
    std::unordered_set<int> ids;
    ids.reserve(static_cast<std::size_t>(p.Num));
    for (int i = 0; i < p.Num; ++i) {
        if (!ids.insert(p.Number[i]).second) {
            throw std::runtime_error(
                "duplicate particle ID " + std::to_string(p.Number[i]) +
                "; CUDA contact history requires stable unique particle IDs");
        }
    }
}

void uploadParticles(const PARTICLE& p, GpuParticleArrays& d)
{
    const int n = p.Num;
    std::vector<double> x(n), y(n), z(n), vx(n), vy(n), vz(n);
    std::vector<double> wx(n), wy(n), wz(n), qw(n), qx(n), qy(n), qz(n);
    std::vector<double> inertiaX(n), inertiaY(n), inertiaZ(n);
    for (int i = 0; i < n; ++i) {
        x[i]=p.Pos[i][0]; y[i]=p.Pos[i][1]; z[i]=p.Pos[i][2];
        vx[i]=p.Vel[i][0]; vy[i]=p.Vel[i][1]; vz[i]=p.Vel[i][2];
        wx[i]=p.AngSpd[i][0]; wy[i]=p.AngSpd[i][1]; wz[i]=p.AngSpd[i][2];
        qw[i]=p.Quat[i][0]; qx[i]=p.Quat[i][1]; qy[i]=p.Quat[i][2]; qz[i]=p.Quat[i][3];
        inertiaX[i]=p.Inertia[i][0]; inertiaY[i]=p.Inertia[i][1];
        inertiaZ[i]=p.Inertia[i][2];
    }
    upload(d.mass, p.Mass, n, "upload particle mass");
    upload(d.inertiaX, inertiaX.data(), n, "upload particle inertia x");
    upload(d.inertiaY, inertiaY.data(), n, "upload particle inertia y");
    upload(d.inertiaZ, inertiaZ.data(), n, "upload particle inertia z");
    upload(d.radius, p.Radius, n, "upload particle radius");
    upload(d.id, p.Number, n, "upload stable particle id");
    upload(d.status, p.Status, n, "upload particle status");
    upload(d.x,x.data(),n,"upload x"); upload(d.y,y.data(),n,"upload y");
    upload(d.z,z.data(),n,"upload z"); upload(d.vx,vx.data(),n,"upload vx");
    upload(d.vy,vy.data(),n,"upload vy"); upload(d.vz,vz.data(),n,"upload vz");
    upload(d.wx,wx.data(),n,"upload wx"); upload(d.wy,wy.data(),n,"upload wy");
    upload(d.wz,wz.data(),n,"upload wz"); upload(d.qw,qw.data(),n,"upload qw");
    upload(d.qx,qx.data(),n,"upload qx"); upload(d.qy,qy.data(),n,"upload qy");
    upload(d.qz,qz.data(),n,"upload qz");
}

double uploadComponents(const PARTICLE& p, GpuComponentArrays& d)
{
    if (p.ComponentNum <= 0 || p.Num <= 0) {
        throw std::runtime_error("CUDA clump topology requires owners and components");
    }
    std::vector<int> order(static_cast<std::size_t>(p.ComponentNum));
    for (int i = 0; i < p.ComponentNum; ++i) order[static_cast<std::size_t>(i)] = i;
    std::sort(order.begin(), order.end(), [&p](int first, int second) {
        if (p.ComponentOwner[first] != p.ComponentOwner[second]) {
            return p.ComponentOwner[first] < p.ComponentOwner[second];
        }
        return p.ComponentId[first] < p.ComponentId[second];
    });
    std::vector<int> id(p.ComponentNum), owner(p.ComponentNum);
    std::vector<int> ownerStart(static_cast<std::size_t>(p.Num) + 1, 0);
    std::vector<double> radius(p.ComponentNum), bodyX(p.ComponentNum),
        bodyY(p.ComponentNum), bodyZ(p.ComponentNum);
    double maximumRadius = 0.0;
    for (int output = 0; output < p.ComponentNum; ++output) {
        const int input = order[static_cast<std::size_t>(output)];
        const int ownerRow = p.ComponentOwner[input];
        if (ownerRow < 0 || ownerRow >= p.Num || p.ComponentId[input] < 0
            || !std::isfinite(p.ComponentRadius[input])
            || p.ComponentRadius[input] <= 0.0) {
            throw std::runtime_error("invalid native clump component topology");
        }
        id[output] = p.ComponentId[input];
        owner[output] = ownerRow;
        radius[output] = p.ComponentRadius[input];
        bodyX[output] = p.ComponentPosBody[input][0];
        bodyY[output] = p.ComponentPosBody[input][1];
        bodyZ[output] = p.ComponentPosBody[input][2];
        ++ownerStart[static_cast<std::size_t>(ownerRow) + 1];
        maximumRadius = std::max(maximumRadius, radius[output]);
    }
    for (int i = 1; i <= p.Num; ++i) {
        ownerStart[static_cast<std::size_t>(i)] +=
            ownerStart[static_cast<std::size_t>(i - 1)];
        const int count = ownerStart[static_cast<std::size_t>(i)]
            - ownerStart[static_cast<std::size_t>(i - 1)];
        if (count < 1 || count > 8) {
            throw std::runtime_error("each CUDA owner requires 1 to 8 components");
        }
    }
    upload(d.id, id.data(), id.size(), "upload component ids");
    upload(d.owner, owner.data(), owner.size(), "upload component owners");
    upload(d.ownerStart, ownerStart.data(), ownerStart.size(),
           "upload component owner starts");
    upload(d.radius, radius.data(), radius.size(), "upload component radii");
    upload(d.bodyX, bodyX.data(), bodyX.size(), "upload component body x");
    upload(d.bodyY, bodyY.data(), bodyY.size(), "upload component body y");
    upload(d.bodyZ, bodyZ.data(), bodyZ.size(), "upload component body z");
    return maximumRadius;
}

void uploadWalls(const WALL& wall, GpuWallArrays& d)
{
    const int n = wall.number;
    std::vector<double> ox(n),oy(n),oz(n),nx(n),ny(n),nz(n);
    for (int i=0;i<n;++i) {
        ox[i]=wall.Orig[i][0]; oy[i]=wall.Orig[i][1]; oz[i]=wall.Orig[i][2];
        nx[i]=wall.N[i][0]; ny[i]=wall.N[i][1]; nz[i]=wall.N[i][2];
    }
    upload(d.ox,ox.data(),n,"upload wall ox"); upload(d.oy,oy.data(),n,"upload wall oy");
    upload(d.oz,oz.data(),n,"upload wall oz"); upload(d.nx,nx.data(),n,"upload wall nx");
    upload(d.ny,ny.data(),n,"upload wall ny"); upload(d.nz,nz.data(),n,"upload wall nz");
}

void uploadBodyStates(const BODYSET& bodies, GpuBodyStateArrays& d)
{
    const int n = bodies.Num;
    std::vector<double> cx(n),cy(n),cz(n),qw(n),qx(n),qy(n),qz(n);
    std::vector<double> vx(n),vy(n),vz(n),wx(n),wy(n),wz(n);
    for (int i=0;i<n;++i) {
        cx[i]=bodies.body[i].MassCenter[0]; cy[i]=bodies.body[i].MassCenter[1];
        cz[i]=bodies.body[i].MassCenter[2];
        qw[i]=bodies.body[i].orien[0]; qx[i]=bodies.body[i].orien[1];
        qy[i]=bodies.body[i].orien[2]; qz[i]=bodies.body[i].orien[3];
        vx[i]=bodies.body[i].Vel[0]; vy[i]=bodies.body[i].Vel[1];
        vz[i]=bodies.body[i].Vel[2];
        wx[i]=bodies.body[i].AngularVel[0]; wy[i]=bodies.body[i].AngularVel[1];
        wz[i]=bodies.body[i].AngularVel[2];
    }
    upload(d.cx,cx.data(),n,"upload body cx"); upload(d.cy,cy.data(),n,"upload body cy");
    upload(d.cz,cz.data(),n,"upload body cz"); upload(d.qw,qw.data(),n,"upload body qw");
    upload(d.qx,qx.data(),n,"upload body qx"); upload(d.qy,qy.data(),n,"upload body qy");
    upload(d.qz,qz.data(),n,"upload body qz"); upload(d.vx,vx.data(),n,"upload body vx");
    upload(d.vy,vy.data(),n,"upload body vy"); upload(d.vz,vz.data(),n,"upload body vz");
    upload(d.wx,wx.data(),n,"upload body wx"); upload(d.wy,wy.data(),n,"upload body wy");
    upload(d.wz,wz.data(),n,"upload body wz");
}

void uploadTriangles(const BODYSET& bodies, GpuTriangleArrays& d)
{
    const int total = countTriangles(bodies);
    std::vector<int> bodyId(total), triId(total);
    std::vector<double> lx1(total),ly1(total),lz1(total),lx2(total),ly2(total),lz2(total);
    std::vector<double> lx3(total),ly3(total),lz3(total),lx4(total),ly4(total),lz4(total);
    std::vector<double> lnx(total),lny(total),lnz(total);
    int output=0;
    for (int body=0;body<bodies.Num;++body) {
        for (int tri=0;tri<bodies.body[body].Size;++tri,++output) {
            const STLTetMesh& t=bodies.body[body].Shape[tri];
            bodyId[output]=body; triId[output]=tri;
            lx1[output]=t.v1[0];ly1[output]=t.v1[1];lz1[output]=t.v1[2];
            lx2[output]=t.v2[0];ly2[output]=t.v2[1];lz2[output]=t.v2[2];
            lx3[output]=t.v3[0];ly3[output]=t.v3[1];lz3[output]=t.v3[2];
            lx4[output]=t.v4[0];ly4[output]=t.v4[1];lz4[output]=t.v4[2];
            lnx[output]=t.normal[0];lny[output]=t.normal[1];lnz[output]=t.normal[2];
        }
    }
    upload(d.bodyId,bodyId.data(),total,"upload triangle body id");
    upload(d.triId,triId.data(),total,"upload triangle local id");
    upload(d.lx1,lx1.data(),total,"upload tri lx1"); upload(d.ly1,ly1.data(),total,"upload tri ly1"); upload(d.lz1,lz1.data(),total,"upload tri lz1");
    upload(d.lx2,lx2.data(),total,"upload tri lx2"); upload(d.ly2,ly2.data(),total,"upload tri ly2"); upload(d.lz2,lz2.data(),total,"upload tri lz2");
    upload(d.lx3,lx3.data(),total,"upload tri lx3"); upload(d.ly3,ly3.data(),total,"upload tri ly3"); upload(d.lz3,lz3.data(),total,"upload tri lz3");
    upload(d.lx4,lx4.data(),total,"upload tri lx4"); upload(d.ly4,ly4.data(),total,"upload tri ly4"); upload(d.lz4,lz4.data(),total,"upload tri lz4");
    upload(d.lnx,lnx.data(),total,"upload tri nx"); upload(d.lny,lny.data(),total,"upload tri ny"); upload(d.lnz,lnz.data(),total,"upload tri nz");
}

} // namespace

class CudaDemSolver::Impl {
public:
    ~Impl()
    {
        gpuFreeTriangleGrid(triangleGrid);
        gpuFreeForces(forces);
        gpuFreeContactHistory(history);
        gpuFreeGrid(particleGrid);
        gpuFreeTriangles(triangles);
        gpuFreeBodyStates(bodies);
        gpuFreeWalls(walls);
        gpuFreeComponents(components);
        gpuFreeParticles(particles);
    }

    GpuParticleArrays particles;
    GpuComponentArrays components;
    GpuWallArrays walls;
    GpuBodyStateArrays bodies;
    GpuTriangleArrays triangles;
    GpuGridArrays particleGrid;
    GpuTriangleGrid triangleGrid;
    GpuContactHistory history;
    GpuForceArrays forces;
    GpuMechanicalParams particleMechanical;
    GpuMechanicalParams bodyMechanical;
    GpuMechanicalParams wallMechanical;
    int stepIndex = 0;
    double meshSize = 0.0;
    double maxRadius = 0.0;
    bool initialized = false;
    std::vector<double> bodyMasses;
    CudaStepStats cachedStats;
    std::size_t freeMemoryBeforeAllocations = 0;
};

CudaDemSolver::CudaDemSolver() : impl_(new Impl) {}
CudaDemSolver::~CudaDemSolver() = default;

void CudaDemSolver::initialize(PARTICLE& p, BODYSET* bodyset,
                               const CONTROL& control, const MECH& particleMech,
                               const MECH& bodyMech, const MECH& wallMech,
                               const ODE& ode)
{
    if (impl_->initialized) throw std::runtime_error("CUDA DEM solver already initialized");
    if (!CudaBackendAvailable()) throw std::runtime_error("no CUDA device is available");
    if (p.Num <= 0) throw std::runtime_error("CUDA DEM solver requires particles");
    validateStableParticleIds(p);
    checkCuda(cudaSetDevice(0), "cudaSetDevice");
    checkCuda(cudaFree(nullptr), "initialize CUDA context");
    std::size_t totalMemory = 0;
    checkCuda(cudaMemGetInfo(&impl_->freeMemoryBeforeAllocations, &totalMemory),
              "query CUDA memory before allocation");

    const int bodyCount = bodyset ? bodyset->Num : 0;
    const int triangleCount = bodyset ? countTriangles(*bodyset) : 0;
    gpuAllocateParticles(impl_->particles, p.Num);
    uploadParticles(p, impl_->particles);
    gpuAllocateComponents(impl_->components, p.ComponentNum, p.Num);
    const double maximumComponentRadius = uploadComponents(p, impl_->components);
    gpuUpdateComponents(impl_->components, impl_->particles);
    gpuAllocateGrid(impl_->particleGrid, p.ComponentNum);
    gpuAllocateContactHistory(impl_->history, p.ComponentNum);
    gpuAllocateForces(impl_->forces, p.Num, bodyCount,
                      control.OutputParticleForce,
                      control.OutputBodyForce);
    if (control.Wall_flag && control.wall.number > 0) {
        gpuAllocateWalls(impl_->walls, control.wall.number);
        uploadWalls(control.wall, impl_->walls);
    }
    if (bodyCount > 0) {
        gpuAllocateBodyStates(impl_->bodies, bodyCount);
        gpuAllocateTriangles(impl_->triangles, triangleCount);
        gpuAllocateTriangleGrid(impl_->triangleGrid, triangleCount);
        uploadBodyStates(*bodyset, impl_->bodies);
        uploadTriangles(*bodyset, impl_->triangles);
        gpuTransformTriangles(impl_->triangles, impl_->bodies);
        impl_->bodyMasses.resize(bodyCount);
        for (int i=0;i<bodyCount;++i) impl_->bodyMasses[i]=bodyset->body[i].Mass;
    }
    impl_->particleMechanical.epsS=particleMech.epsS;
    impl_->particleMechanical.mu=particleMech.mu;
    impl_->particleMechanical.epsN=particleMech.epsN;
    impl_->particleMechanical.cohesion=particleMech.c;
    impl_->particleMechanical.beta=particleMech.beta;
    impl_->particleMechanical.muT=particleMech.mu_T;
    impl_->particleMechanical.muR=particleMech.mu_R;
    impl_->particleMechanical.kN=particleMech.kN;
    impl_->particleMechanical.dt=ode.StepSize;
    impl_->particleMechanical.gx=control.g[0];
    impl_->particleMechanical.gy=control.g[1];
    impl_->particleMechanical.gz=control.g[2];
    impl_->particleMechanical.gravityEnabled=control.Grav_flag ? 1 : 0;
    impl_->particleMechanical.universalGravityEnabled=
        control.Universal_gravitation_flag ? 1 : 0;
    impl_->bodyMechanical.epsS=bodyMech.epsS;
    impl_->bodyMechanical.mu=bodyMech.mu;
    impl_->bodyMechanical.epsN=bodyMech.epsN;
    impl_->bodyMechanical.cohesion=bodyMech.c;
    impl_->bodyMechanical.beta=bodyMech.beta;
    impl_->bodyMechanical.muT=bodyMech.mu_T;
    impl_->bodyMechanical.muR=bodyMech.mu_R;
    impl_->bodyMechanical.kN=bodyMech.kN;
    impl_->bodyMechanical.dt=ode.StepSize;
    impl_->wallMechanical.epsS=wallMech.epsS;
    impl_->wallMechanical.mu=wallMech.mu;
    impl_->wallMechanical.epsN=wallMech.epsN;
    impl_->wallMechanical.cohesion=wallMech.c;
    impl_->wallMechanical.beta=wallMech.beta;
    impl_->wallMechanical.muT=wallMech.mu_T;
    impl_->wallMechanical.muR=wallMech.mu_R;
    impl_->wallMechanical.kN=wallMech.kN;
    impl_->wallMechanical.dt=ode.StepSize;
    impl_->meshSize=2.5 * maximumComponentRadius;
    impl_->maxRadius=maximumComponentRadius;
    impl_->initialized=true;
}

void CudaDemSolver::computeForces(double, BODYSET* bodyset, double** bodyForces,
                                  bool collectDiagnostics)
{
    if (!impl_->initialized) throw std::runtime_error("CUDA DEM solver is not initialized");
    if (bodyset && bodyset->Num > 0) {
        uploadBodyStates(*bodyset, impl_->bodies);
        gpuTransformTriangles(impl_->triangles, impl_->bodies);
    }
    gpuUpdateComponents(impl_->components, impl_->particles);
    gpuBuildComponentGrid(impl_->particleGrid, impl_->components, impl_->meshSize);
    if (impl_->triangles.n > 0) {
        gpuBuildTriangleGrid(impl_->triangleGrid, impl_->triangles,
                             impl_->meshSize, impl_->maxRadius);
    }
    gpuClearForces(impl_->forces);
    const int contactStep = impl_->stepIndex++;
    GpuMechanicalParams particleMechanical = impl_->particleMechanical;
    GpuMechanicalParams bodyMechanical = impl_->bodyMechanical;
    GpuMechanicalParams wallMechanical = impl_->wallMechanical;
    // The pre-integration force evaluation initializes the leapfrog state but
    // advances no physical time.  Keep dashpot forces active while preventing
    // an unphysical full-dt increment of tangent/roll/twist history.
    if (contactStep == 0) {
        particleMechanical.dt = 0.0;
        bodyMechanical.dt = 0.0;
        wallMechanical.dt = 0.0;
    }
    gpuComputeContacts(impl_->particles, impl_->components,
                       impl_->particleGrid, impl_->walls,
                       impl_->bodies, impl_->triangles, impl_->triangleGrid,
                       impl_->history, impl_->forces,
                       particleMechanical, bodyMechanical, wallMechanical,
                       contactStep);
    gpuComputeParticleAcceleration(impl_->particles);

    if (collectDiagnostics) {
    CudaStepStats stats;
    download(&stats.particleContacts, impl_->forces.particleContactCount, 1,
             "download particle contact count");
    download(&stats.bodyContacts, impl_->forces.bodyContactCount, 1,
             "download body contact count");
    download(&stats.bodyMaxDepth, impl_->forces.bodyMaxDepth, 1,
             "download body max depth");
    download(&stats.bodyContactFz, impl_->forces.bodyContactFz, 1,
             "download body contact Fz");
    download(&stats.particleHistoryOverflows,
             impl_->forces.particleHistoryOverflowCount, 1,
             "download particle history overflow count");
    download(&stats.bodyHistoryOverflows,
             impl_->forces.bodyHistoryOverflowCount, 1,
             "download body history overflow count");
    download(&stats.particleHistoryHighWater,
             impl_->forces.particleHistoryHighWater, 1,
             "download particle history high-water mark");
    download(&stats.bodyHistoryHighWater,
             impl_->forces.bodyHistoryHighWater, 1,
             "download body history high-water mark");
    if (impl_->triangleGrid.overflowCount) {
        download(&stats.triangleGridOverflows, impl_->triangleGrid.overflowCount, 1,
                 "download triangle grid overflow count");
    }
    impl_->cachedStats=stats;
    if (stats.particleHistoryOverflows != 0) {
        throw std::runtime_error("CUDA particle contact-history capacity exceeded; increase kGpuParticleHistorySlots");
    }
    if (stats.bodyHistoryOverflows != 0) {
        throw std::runtime_error("CUDA wall/body contact-history capacity exceeded; increase the configured history slots");
    }
    if (stats.triangleGridOverflows != 0) {
        throw std::runtime_error("CUDA triangle-grid coordinate/count overflow");
    }
    }

    if (bodyForces && impl_->bodies.n > 0) {
        std::vector<double> result(static_cast<std::size_t>(impl_->bodies.n)*6);
        download(result.data(), impl_->forces.bodyResult, result.size(),
                 "download body resultant forces");
        for (int body=0;body<impl_->bodies.n;++body) {
            for (int k=0;k<6;++k) bodyForces[body][k]=result[body*6+k];
            if (impl_->particleMechanical.gravityEnabled) {
                bodyForces[body][0] += impl_->particleMechanical.gx*impl_->bodyMasses[body];
                bodyForces[body][1] += impl_->particleMechanical.gy*impl_->bodyMasses[body];
                bodyForces[body][2] += impl_->particleMechanical.gz*impl_->bodyMasses[body];
            }
        }
    }
}

void CudaDemSolver::initializeParticleHalfStep(double halfDt)
{
    gpuInitializeParticleLeapfrog(impl_->particles, halfDt);
}

double CudaDemSolver::advanceParticles(double fullDt)
{
    gpuAdvanceParticlePosition(impl_->particles, fullDt);
    // The CUDA broad phase is rebuilt every step, so the legacy CPU tree
    // displacement threshold is unnecessary and would force a device-to-host
    // reduction/synchronization in the hot loop.
    return 0.0;
}

void CudaDemSolver::finishParticleStep(double fullDt, double)
{
    gpuAdvanceParticleHalfVelocity(impl_->particles, fullDt);
}

void CudaDemSolver::saveContactHistory(const std::string& filename) const
{
    if (!impl_->initialized) {
        throw std::runtime_error("cannot save contact history before CUDA initialization");
    }
    const std::size_t particleSlots = static_cast<std::size_t>(impl_->history.particleCount)
        * kGpuParticleHistorySlots;
    const std::size_t wallSlots = static_cast<std::size_t>(impl_->history.particleCount)
        * kGpuWallHistorySlots;
    const std::size_t bodySlots = static_cast<std::size_t>(impl_->history.particleCount)
        * kGpuBodyHistorySlots;
    std::vector<std::uint64_t> particleKeys(particleSlots), wallKeys(wallSlots),
        bodyKeys(bodySlots);
    std::vector<double> particleValues(particleSlots * 9), wallValues(wallSlots * 9),
        bodyValues(bodySlots * 9);
    download(particleKeys.data(), impl_->history.particleKeys, particleSlots,
             "download particle history keys");
    download(wallKeys.data(), impl_->history.wallKeys, wallSlots,
             "download wall history keys");
    download(bodyKeys.data(), impl_->history.bodyKeys, bodySlots,
             "download body history keys");
    download(particleValues.data(), impl_->history.particleValues, particleValues.size(),
             "download particle history values");
    download(wallValues.data(), impl_->history.wallValues, wallValues.size(),
             "download wall history values");
    download(bodyValues.data(), impl_->history.bodyValues, bodyValues.size(),
             "download body history values");

    std::ofstream output(filename, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("cannot create CUDA contact-history checkpoint: " + filename);
    const std::array<char, 8> magic{{'D','M','H','I','S','T','1','\0'}};
    const std::uint32_t version = 1;
    const std::uint32_t particleCount = static_cast<std::uint32_t>(impl_->history.particleCount);
    const std::array<std::uint32_t, 3> slotCounts{{
        kGpuParticleHistorySlots, kGpuWallHistorySlots, kGpuBodyHistorySlots}};
    const auto write = [&output](const auto* values, std::size_t count) {
        output.write(reinterpret_cast<const char*>(values),
                     static_cast<std::streamsize>(sizeof(*values) * count));
    };
    write(magic.data(), magic.size());
    write(&version, 1);
    write(&particleCount, 1);
    write(slotCounts.data(), slotCounts.size());
    write(particleKeys.data(), particleKeys.size());
    write(wallKeys.data(), wallKeys.size());
    write(bodyKeys.data(), bodyKeys.size());
    write(particleValues.data(), particleValues.size());
    write(wallValues.data(), wallValues.size());
    write(bodyValues.data(), bodyValues.size());
    if (!output) throw std::runtime_error("failed to write CUDA contact-history checkpoint: " + filename);
}

void CudaDemSolver::loadContactHistory(const std::string& filename)
{
    if (!impl_->initialized) {
        throw std::runtime_error("cannot load contact history before CUDA initialization");
    }
    std::ifstream input(filename, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open CUDA contact-history checkpoint: " + filename);
    const auto read = [&input, &filename](auto* values, std::size_t count) {
        input.read(reinterpret_cast<char*>(values),
                   static_cast<std::streamsize>(sizeof(*values) * count));
        if (!input) throw std::runtime_error("truncated CUDA contact-history checkpoint: " + filename);
    };
    std::array<char, 8> magic{};
    std::uint32_t version = 0;
    std::uint32_t particleCount = 0;
    std::array<std::uint32_t, 3> slotCounts{};
    read(magic.data(), magic.size());
    read(&version, 1);
    read(&particleCount, 1);
    read(slotCounts.data(), slotCounts.size());
    const std::array<char, 8> expectedMagic{{'D','M','H','I','S','T','1','\0'}};
    const std::array<std::uint32_t, 3> expectedSlots{{
        kGpuParticleHistorySlots, kGpuWallHistorySlots, kGpuBodyHistorySlots}};
    if (magic != expectedMagic || version != 1
        || particleCount != static_cast<std::uint32_t>(impl_->history.particleCount)
        || slotCounts != expectedSlots) {
        throw std::runtime_error("incompatible CUDA contact-history checkpoint: " + filename);
    }
    const std::size_t particleSlots = static_cast<std::size_t>(particleCount)
        * kGpuParticleHistorySlots;
    const std::size_t wallSlots = static_cast<std::size_t>(particleCount)
        * kGpuWallHistorySlots;
    const std::size_t bodySlots = static_cast<std::size_t>(particleCount)
        * kGpuBodyHistorySlots;
    std::vector<std::uint64_t> particleKeys(particleSlots), wallKeys(wallSlots),
        bodyKeys(bodySlots);
    std::vector<double> particleValues(particleSlots * 9), wallValues(wallSlots * 9),
        bodyValues(bodySlots * 9);
    read(particleKeys.data(), particleKeys.size());
    read(wallKeys.data(), wallKeys.size());
    read(bodyKeys.data(), bodyKeys.size());
    read(particleValues.data(), particleValues.size());
    read(wallValues.data(), wallValues.size());
    read(bodyValues.data(), bodyValues.size());
    char trailing = 0;
    if (input.read(&trailing, 1)) {
        throw std::runtime_error("oversized CUDA contact-history checkpoint: " + filename);
    }
    upload(impl_->history.particleKeys, particleKeys.data(), particleKeys.size(),
           "upload particle history keys");
    upload(impl_->history.wallKeys, wallKeys.data(), wallKeys.size(),
           "upload wall history keys");
    upload(impl_->history.bodyKeys, bodyKeys.data(), bodyKeys.size(),
           "upload body history keys");
    upload(impl_->history.particleValues, particleValues.data(), particleValues.size(),
           "upload particle history values");
    upload(impl_->history.wallValues, wallValues.data(), wallValues.size(),
           "upload wall history values");
    upload(impl_->history.bodyValues, bodyValues.data(), bodyValues.size(),
           "upload body history values");
    checkCuda(cudaMemset(impl_->history.particleLastSeen, 0xff,
                         sizeof(int) * particleSlots), "reset particle history age");
    checkCuda(cudaMemset(impl_->history.wallLastSeen, 0xff,
                         sizeof(int) * wallSlots), "reset wall history age");
    checkCuda(cudaMemset(impl_->history.bodyLastSeen, 0xff,
                         sizeof(int) * bodySlots), "reset body history age");
    impl_->stepIndex = 0;
}

void CudaDemSolver::downloadParticleState(PARTICLE& p) const
{
    const int n=impl_->particles.n;
    if (p.Num != n) throw std::runtime_error("particle count changed during CUDA simulation");
    std::vector<double> x(n),y(n),z(n),vx(n),vy(n),vz(n),wx(n),wy(n),wz(n);
    std::vector<double> qw(n),qx(n),qy(n),qz(n);
    download(x.data(),impl_->particles.x,n,"download x"); download(y.data(),impl_->particles.y,n,"download y");
    download(z.data(),impl_->particles.z,n,"download z"); download(vx.data(),impl_->particles.vx,n,"download vx");
    download(vy.data(),impl_->particles.vy,n,"download vy"); download(vz.data(),impl_->particles.vz,n,"download vz");
    download(wx.data(),impl_->particles.wx,n,"download wx"); download(wy.data(),impl_->particles.wy,n,"download wy");
    download(wz.data(),impl_->particles.wz,n,"download wz"); download(qw.data(),impl_->particles.qw,n,"download qw");
    download(qx.data(),impl_->particles.qx,n,"download qx"); download(qy.data(),impl_->particles.qy,n,"download qy");
    download(qz.data(),impl_->particles.qz,n,"download qz");
    for(int i=0;i<n;++i) {
        p.Pos[i]=vector3d(x[i],y[i],z[i]); p.Vel[i]=vector3d(vx[i],vy[i],vz[i]);
        p.AngSpd[i]=vector3d(wx[i],wy[i],wz[i]); p.Quat[i]=vector4d(qw[i],qx[i],qy[i],qz[i]);
    }
}

void CudaDemSolver::downloadForceOutput(FORCE& f, bool particleForces,
                                        bool bodyForces) const
{
    const int n=impl_->forces.particleCount;
    const int nb=impl_->forces.bodyCount;
    const std::size_t p6=static_cast<std::size_t>(n)*6;
    if (particleForces) {
        std::vector<double> contact(p6),field(p6),wall(p6),total(p6);
        download(contact.data(),impl_->forces.contact,p6,"download contact force output");
        download(field.data(),impl_->forces.field,p6,"download field force output");
        download(wall.data(),impl_->forces.wall,p6,"download wall force output");
        download(total.data(),impl_->forces.total,p6,"download total force output");
        for(int p=0;p<n;++p) for(int k=0;k<6;++k) {
            f.Fc[p][k]=contact[p*6+k]; f.Ff[p][k]=field[p*6+k];
            f.Fw[p][k]=wall[p*6+k]; f.Fd[p][k]=0.0; f.FR[p][k]=total[p*6+k];
        }
    }
    if(nb>0 && (particleForces || bodyForces)) {
        std::vector<double> bp;
        std::vector<double> detail;
        if (particleForces) {
            bp.resize(static_cast<std::size_t>(nb)*p6);
            download(bp.data(),impl_->forces.bodyParticle,bp.size(),
                     "download body-particle force output");
        }
        if (bodyForces) {
            detail.resize(static_cast<std::size_t>(nb)*n*19);
            download(detail.data(),impl_->forces.bodyDetail,detail.size(),
                     "download body detail output");
        }
        for(int b=0;b<nb;++b) for(int p=0;p<n;++p) {
            const int bpBase=(b*n+p)*6;
            if (particleForces) {
                for(int k=0;k<6;++k) f.FI_pt[b][p][k]=bp[bpBase+k];
            }
            if (bodyForces) {
                const int detailBase=(b*n+p)*19;
                for(int k=0;k<19;++k) f.ForceBody[b][p][k]=detail[detailBase+k];
            }
        }
    }
    f.bodyContactCount=impl_->cachedStats.bodyContacts;
    f.bodyContactMaxDepth=impl_->cachedStats.bodyMaxDepth;
    f.bodyContactFz=impl_->cachedStats.bodyContactFz;
}

CudaStepStats CudaDemSolver::stepStats() const { return impl_->cachedStats; }
std::size_t CudaDemSolver::particleGridEntryCount() const { return impl_->particleGrid.n; }
std::size_t CudaDemSolver::triangleGridEntryCount() const { return impl_->triangleGrid.entryCount; }
std::size_t CudaDemSolver::deviceMemoryBytes() const
{
    std::size_t freeMemory = 0;
    std::size_t totalMemory = 0;
    checkCuda(cudaMemGetInfo(&freeMemory, &totalMemory), "query CUDA memory usage");
    return impl_->freeMemoryBeforeAllocations > freeMemory
        ? impl_->freeMemoryBeforeAllocations - freeMemory : 0;
}

bool CudaBackendAvailable()
{
    int count=0;
    return cudaGetDeviceCount(&count)==cudaSuccess && count>0;
}

const char* CudaDeviceDescription()
{
    static std::array<char,256> description{};
    if (description[0]=='\0') {
        cudaDeviceProp property{};
        if (cudaGetDeviceProperties(&property,0)==cudaSuccess) {
            std::snprintf(description.data(),description.size(),
                          "%s, compute capability %d.%d, %.1f GiB",
                          property.name,property.major,property.minor,
                          static_cast<double>(property.totalGlobalMem)/(1024.0*1024.0*1024.0));
        } else {
            std::snprintf(description.data(),description.size(),"unavailable");
        }
    }
    return description.data();
}
