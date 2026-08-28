#pragma once

#include "force.h"
#include "multibody.h"
#include "particle.h"

#include <array>
#include <memory>
#include <vector>

struct AsyncOutputOptions {
    bool particleState = false;
    bool particleForce = false;
    bool bodyState = false;
    bool bodyForce = false;
    bool any(int bodyCount) const {
        return particleState || particleForce ||
               (bodyCount > 0 && (bodyState || bodyForce));
    }
};

class OutputFrame {
public:
    OutputFrame(const PARTICLE&, int, const AsyncOutputOptions&);
    ~OutputFrame();
    OutputFrame(const OutputFrame&) = delete;
    OutputFrame& operator=(const OutputFrame&) = delete;
    bool needsParticleState() const { return options_.particleState; }
    bool needsParticleForce() const { return options_.particleForce; }
    bool needsBodyForce() const { return options_.bodyForce && bodyCount_ > 0; }
    PARTICLE& particles();
    FORCE& forces();
    void captureBodies(const BODYSET&, double**, double**);
    void setMetadata(int, double);

private:
    friend class AsyncOutputWriter;
    void writeFiles();
    AsyncOutputOptions options_;
    std::unique_ptr<PARTICLE> particles_;
    std::unique_ptr<FORCE> forces_;
    BODYSET bodies_;
    std::unique_ptr<BODY[]> bodyStorage_;
    std::vector<std::array<double, 6>> bodyForces_;
    std::vector<double*> bodyForceRows_;
    std::vector<std::array<double, 6>> bodyImpulses_;
    std::vector<double*> bodyImpulseRows_;
    int bodyCount_ = 0;
    int step_ = -1;
    double time_ = 0.0;
};

class AsyncOutputWriter {
public:
    AsyncOutputWriter(const PARTICLE&, int, const AsyncOutputOptions&);
    ~AsyncOutputWriter();
    AsyncOutputWriter(const AsyncOutputWriter&) = delete;
    AsyncOutputWriter& operator=(const AsyncOutputWriter&) = delete;
    OutputFrame& acquire();
    void submit(OutputFrame&);
    void finish();
private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};
