#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

class PARTICLE;
class BODYSET;
class FORCE;
struct CONTROL;
struct MECH;
struct ODE;

struct CudaStepStats {
    int particleContacts = 0;
    int bodyContacts = 0;
    double bodyMaxDepth = 0.0;
    double bodyContactFz = 0.0;
    int particleHistoryOverflows = 0;
    int bodyHistoryOverflows = 0;
    int particleHistoryHighWater = 0;
    int bodyHistoryHighWater = 0;
    int triangleGridOverflows = 0;
    // Per-body PT component-manifold counts.  Columns are active, primary
    // stick, primary slide, primary twist, slide flag, twist flag, and both
    // flags.  Primary classes use the fixed priority twist > slide > stick.
    std::vector<std::array<std::uint64_t, 7>> bodyRegimeInstant;
    std::vector<std::array<std::uint64_t, 7>> bodyRegimeInterval;
    std::uint64_t bodyRegimeIntervalNativeSteps = 0;
};

// Persistent CUDA DEM solver. Particle state, contact history, spatial grids,
// and force buffers remain device-resident for the complete simulation.
class CudaDemSolver {
public:
    CudaDemSolver();
    ~CudaDemSolver();
    CudaDemSolver(const CudaDemSolver&) = delete;
    CudaDemSolver& operator=(const CudaDemSolver&) = delete;

    void initialize(PARTICLE& particles, BODYSET* bodies,
                    const CONTROL& control, const MECH& particleMech,
                    const MECH& bodyMech, const MECH& wallMech,
                    const ODE& ode);
    void computeForces(double time, BODYSET* bodies, double** bodyForces,
                       bool collectDiagnostics);
    void initializeParticleHalfStep(double halfDt);
    double advanceParticles(double fullDt);
    void finishParticleStep(double fullDt, double halfDt);
    void loadContactHistory(const std::string& filename);
    void saveContactHistory(const std::string& filename) const;
    bool hasRestoredParticleLeapfrogState() const;

    void downloadParticleState(PARTICLE& particles) const;
    void downloadForceOutput(FORCE& force, bool particleForces = true,
                             bool bodyForces = true) const;
    CudaStepStats stepStats() const;
    std::size_t particleGridEntryCount() const;
    std::size_t triangleGridEntryCount() const;
    std::size_t deviceMemoryBytes() const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

bool CudaBackendAvailable();
const char* CudaDeviceDescription();
