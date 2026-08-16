#include "simulation.h"

#include "async_output.h"
#include "cuda_backend.h"
#include "dyneq.h"

#include <algorithm>
#include <array>
#include <iostream>
#include <memory>
#include <vector>

namespace {

class BodyForceStorage {
public:
    explicit BodyForceStorage(int count) : rows_(count), pointers_(count) {
        for (int i = 0; i < count; ++i) pointers_[i] = rows_[i].data();
    }
    double** data() { return pointers_.data(); }

private:
    std::vector<std::array<double, 6>> rows_;
    std::vector<double*> pointers_;
};

class BodyAccelerationStorage {
public:
    explicit BodyAccelerationStorage(int count) : rows_(count), pointers_(count) {
        for (int i = 0; i < count; ++i) pointers_[i] = rows_[i].data();
    }
    double** data() { return pointers_.data(); }

private:
    std::vector<std::array<double, 6>> rows_;
    std::vector<double*> pointers_;
};

AsyncOutputOptions outputOptions(const CONTROL& control) {
    return {control.OutputParticleState, control.OutputParticleForce,
            control.OutputBodyState, control.OutputBodyForce};
}

void enqueueCudaFrame(AsyncOutputWriter* writer, CudaDemSolver& solver,
                      BODYSET* bodies, double** bodyForces,
                      int step, double time) {
    if (!writer) return;
    OutputFrame& frame = writer->acquire();
    if (frame.needsParticleState()) solver.downloadParticleState(frame.particles());
    if (frame.needsParticleForce() || frame.needsBodyForce()) {
        solver.downloadForceOutput(frame.forces(), frame.needsParticleForce(),
                                   frame.needsBodyForce());
    }
    if (bodies) frame.captureBodies(*bodies, bodyForces);
    frame.setMetadata(step, time);
    writer->submit(frame);
}

}  // namespace

void Simulation::IntegrateDem() {
    const double halfDt = 0.5 * ode.StepSize;
    CudaDemSolver solver;
    solver.initialize(pt, nullptr, control, mechPP, mechPT, mechPW, ode);
    const AsyncOutputOptions options = outputOptions(control);
    std::unique_ptr<AsyncOutputWriter> output;
    if (options.any(0)) output = std::make_unique<AsyncOutputWriter>(pt, 0, options);

    std::cout << "CUDA DEM solver: " << CudaDeviceDescription() << std::endl;
    solver.computeForces(0.0, nullptr, nullptr, true);
    std::cout << "Persistent CUDA memory = "
              << solver.deviceMemoryBytes() / (1024.0 * 1024.0) << " MiB" << std::endl;
    enqueueCudaFrame(output.get(), solver, nullptr, nullptr, -1, 0.0);
    solver.initializeParticleHalfStep(halfDt);

    int outputCounter = 0;
    for (int step = ode.StartStep; step < ode.EndStep; ++step) {
        const double nextTime = (step + 1) * ode.StepSize;
        pt.Vmax = solver.advanceParticles(ode.StepSize);
        const bool outputDue = outputCounter + 1 == ode.OutputInterval ||
                               step == ode.EndStep - 1;
        solver.computeForces(nextTime, nullptr, nullptr, outputDue);
        solver.finishParticleStep(ode.StepSize, halfDt);

        ++outputCounter;
        if (outputCounter == ode.OutputInterval || step == ode.EndStep - 1) {
            enqueueCudaFrame(output.get(), solver, nullptr, nullptr,
                             step, nextTime);
            const auto stats = solver.stepStats();
            std::cout << "Step = " << step + 1 << "/" << ode.EndStep
                      << ", particle grid entries = "
                      << solver.particleGridEntryCount()
                      << ", contacts = " << stats.particleContacts << std::endl;
            outputCounter = 0;
        }
    }
    if (output) output->finish();
}

void Simulation::IntegrateDemMultiBody() {
    bodyset.updateMatrix();
    const double halfDt = 0.5 * ode.StepSize;
    BodyForceStorage bodyForceStorage(bodyset.Num);
    BodyAccelerationStorage bodyAccelerationStorage(bodyset.Num);
    double** bodyForces = bodyForceStorage.data();
    double** bodyAccelerations = bodyAccelerationStorage.data();
    VHALF bodyHalf(bodyset.Num);
    ACC bodyAcceleration(bodyset.Num);

    CudaDemSolver solver;
    solver.initialize(pt, &bodyset, control, mechPP, mechPT, mechPW, ode);
    const AsyncOutputOptions options = outputOptions(control);
    std::unique_ptr<AsyncOutputWriter> output;
    if (options.any(bodyset.Num))
        output = std::make_unique<AsyncOutputWriter>(pt, bodyset.Num, options);

    std::cout << "CUDA DEM-MBD solver: " << CudaDeviceDescription() << std::endl;
    solver.computeForces(0.0, &bodyset, bodyForces, true);
    std::cout << "Persistent CUDA memory = "
              << solver.deviceMemoryBytes() / (1024.0 * 1024.0) << " MiB" << std::endl;
    enqueueCudaFrame(output.get(), solver, &bodyset, bodyForces, -1, 0.0);
    solver.initializeParticleHalfStep(halfDt);

    DynEqnBodySet(0.0, bodyset, bodyForces, bodyAccelerations);
    for (int body = 0; body < bodyset.Num; ++body) {
        const int state = bodyset.body[body].state;
        for (int axis = 0; axis < 3; ++axis) {
            bodyAcceleration.acc[body][axis] = bodyAccelerations[body][axis];
            bodyAcceleration.domg[body][axis] = bodyAccelerations[body][axis + 3];
        }
        if (state == BODY_DYNAMIC) {
            bodyHalf.Vel[body] = bodyset.body[body].Vel +
                                 bodyAcceleration.acc[body] * halfDt;
            bodyHalf.AngSpd[body] = bodyset.body[body].AngularVel +
                                    bodyAcceleration.domg[body] * halfDt;
        } else if (state == BODY_FIXED) {
            bodyAcceleration.acc[body] = vector3d{0.0, 0.0, 0.0};
            bodyAcceleration.domg[body] = vector3d{0.0, 0.0, 0.0};
            bodyHalf.Vel[body] = vector3d{0.0, 0.0, 0.0};
            bodyHalf.AngSpd[body] = vector3d{0.0, 0.0, 0.0};
        } else {
            bodyAcceleration.acc[body] = vector3d{0.0, 0.0, 0.0};
            bodyAcceleration.domg[body] = vector3d{0.0, 0.0, 0.0};
            bodyHalf.Vel[body] = bodyset.body[body].Vel;
            bodyHalf.AngSpd[body] = bodyset.body[body].AngularVel;
        }
        bodyHalf.quat[body] = (bodyset.body[body].orien +
            gama(bodyset.body[body].orien, bodyset.body[body].AngularVel) * halfDt).normalize();
        bodyHalf.dQuat[body] = gama(bodyHalf.quat[body], bodyHalf.AngSpd[body]);
    }

    int outputCounter = 0;
    for (int step = ode.StartStep; step < ode.EndStep; ++step) {
        const double currentTime = step * ode.StepSize;
        const double nextTime = (step + 1) * ode.StepSize;
        pt.Vmax = solver.advanceParticles(ode.StepSize);

        for (int body = 0; body < bodyset.Num; ++body) {
            if (bodyset.body[body].state == BODY_FIXED) continue;
            bodyset.body[body].MassCenter += bodyHalf.Vel[body] * ode.StepSize;
            bodyset.body[body].orien = (bodyset.body[body].orien +
                bodyHalf.dQuat[body] * ode.StepSize).normalize();
            bodyset.body[body].Vel = bodyHalf.Vel[body] +
                                     bodyAcceleration.acc[body] * halfDt;
            bodyset.body[body].AngularVel = bodyHalf.AngSpd[body] +
                                            bodyAcceleration.domg[body] * halfDt;
        }
        bodyset.updateMatrix();

        const bool outputDue = outputCounter + 1 == ode.OutputInterval ||
                               step == ode.EndStep - 1;
        solver.computeForces(nextTime, &bodyset, bodyForces, outputDue);
        solver.finishParticleStep(ode.StepSize, halfDt);
        DynEqnBodySet(nextTime, bodyset, bodyForces, bodyAccelerations);
        for (int body = 0; body < bodyset.Num; ++body) {
            const int state = bodyset.body[body].state;
            for (int axis = 0; axis < 3; ++axis) {
                bodyAcceleration.acc[body][axis] = bodyAccelerations[body][axis];
                bodyAcceleration.domg[body][axis] = bodyAccelerations[body][axis + 3];
            }
            if (state == BODY_DYNAMIC) {
                bodyHalf.Vel[body] += bodyAcceleration.acc[body] * ode.StepSize;
                bodyHalf.AngSpd[body] += bodyAcceleration.domg[body] * ode.StepSize;
            } else if (state == BODY_FIXED) {
                bodyset.body[body].Vel = vector3d{0.0, 0.0, 0.0};
                bodyset.body[body].AngularVel = vector3d{0.0, 0.0, 0.0};
                bodyHalf.Vel[body] = vector3d{0.0, 0.0, 0.0};
                bodyHalf.AngSpd[body] = vector3d{0.0, 0.0, 0.0};
            } else {
                // Prescribed motion currently means constant linear/angular
                // velocity from the CSV; computed contact loads do not alter it.
                bodyHalf.Vel[body] = bodyset.body[body].Vel;
                bodyHalf.AngSpd[body] = bodyset.body[body].AngularVel;
            }
            bodyHalf.quat[body] = (bodyset.body[body].orien +
                gama(bodyset.body[body].orien, bodyset.body[body].AngularVel) * halfDt).normalize();
            bodyHalf.dQuat[body] = gama(bodyHalf.quat[body], bodyHalf.AngSpd[body]);
        }

        ++outputCounter;
        if (outputCounter == ode.OutputInterval || step == ode.EndStep - 1) {
            enqueueCudaFrame(output.get(), solver, &bodyset, bodyForces,
                             step, currentTime);
            const auto stats = solver.stepStats();
            std::cout << "  " << step + 1 << "/" << ode.EndStep
                      << ", particle/triangle grid = "
                      << solver.particleGridEntryCount() << "/"
                      << solver.triangleGridEntryCount()
                      << ", body contacts = " << stats.bodyContacts << std::endl;
            outputCounter = 0;
        }
    }
    if (output) output->finish();
}
