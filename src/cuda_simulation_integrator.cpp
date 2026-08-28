#include "simulation.h"

#include "async_output.h"
#include "cuda_backend.h"
#include "dyneq.h"
#include "file_utils.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void restoreContactHistoryIfPresent(CudaDemSolver& solver) {
    const auto input = file_utils::resolveCasePath("InputFile/contact_history.bt");
    std::error_code error;
    const bool exists = std::filesystem::exists(input, error);
    if (error) {
        throw file_utils::pathError("Cannot inspect / 无法检查",
            "CUDA contact-history checkpoint", input, error.message());
    }
    if (exists) solver.loadContactHistory(input.string());
}

void saveContactHistory(CudaDemSolver& solver) {
    const auto output = file_utils::prepareOutputFile(
        "OutputFile/contact_history.bt", "CUDA contact-history checkpoint");
    solver.saveContactHistory(output.string());
}

class BodyForceStorage {
public:
    explicit BodyForceStorage(int count) : rows_(count), pointers_(count) {
        for (int i = 0; i < count; ++i) pointers_[i] = rows_[i].data();
    }
    double** data() { return pointers_.data(); }
    void clear() {
        for (auto& row : rows_) row.fill(0.0);
    }
    void add(double** source, double scale) {
        for (std::size_t body = 0; body < rows_.size(); ++body) {
            for (int axis = 0; axis < 6; ++axis) {
                rows_[body][axis] += source[body][axis] * scale;
            }
        }
    }

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
                      BODYSET* bodies, double** bodyForces, double** bodyImpulses,
                      int step, double time) {
    if (!writer) return;
    OutputFrame& frame = writer->acquire();
    if (frame.needsParticleState()) solver.downloadParticleState(frame.particles());
    if (frame.needsParticleForce() || frame.needsBodyForce()) {
        solver.downloadForceOutput(frame.forces(), frame.needsParticleForce(),
                                   frame.needsBodyForce());
    }
    if (bodies) frame.captureBodies(*bodies, bodyForces, bodyImpulses);
    frame.setMetadata(step, time);
    writer->submit(frame);
}

struct PrescribedMotionKey {
    int step = 0;
    vector3d linear{0.0, 0.0, 0.0};
    vector3d angular{0.0, 0.0, 0.0};
};

std::vector<PrescribedMotionKey> loadPrescribedMotion(
    const std::string& filename, int startStep, int endStep)
{
    if (filename.empty()) return {};
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("cannot open prescribed motion file: " + filename);
    }
    std::vector<PrescribedMotionKey> keys;
    std::string line;
    int lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        const auto comment = line.find('#');
        if (comment != std::string::npos) line.erase(comment);
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        PrescribedMotionKey key;
        if (!(row >> key.step)) continue;
        if (!(row >> key.linear[0] >> key.linear[1] >> key.linear[2]
                  >> key.angular[0] >> key.angular[1] >> key.angular[2])) {
            throw std::runtime_error(
                "invalid prescribed motion row " + std::to_string(lineNumber));
        }
        std::string trailing;
        if (row >> trailing) {
            throw std::runtime_error(
                "extra prescribed motion column at row " + std::to_string(lineNumber));
        }
        if ((!keys.empty() && key.step <= keys.back().step)
            || key.step < startStep || key.step >= endStep) {
            throw std::runtime_error(
                "prescribed motion steps must increase within the integration range");
        }
        keys.push_back(key);
    }
    if (keys.empty() || keys.front().step != startStep) {
        throw std::runtime_error(
            "prescribed motion must begin at ODE_StartStep");
    }
    return keys;
}

void applyPrescribedMotion(
    const std::vector<PrescribedMotionKey>& keys, std::size_t& cursor,
    int step, double halfDt, BODYSET& bodies, VHALF* half)
{
    if (keys.empty()) return;
    while (cursor + 1 < keys.size() && keys[cursor + 1].step <= step) ++cursor;
    if (bodies.Num != 1 || bodies.body[0].state != BODY_PRESCRIBED) {
        throw std::runtime_error(
            "prescribed motion schedule currently requires one prescribed body");
    }
    BODY& body = bodies.body[0];
    body.Vel = keys[cursor].linear;
    body.AngularVel = keys[cursor].angular;
    if (half) {
        half->Vel[0] = body.Vel;
        half->AngSpd[0] = body.AngularVel;
        half->quat[0] = (body.orien
            + gama(body.orien, body.AngularVel) * halfDt).normalize();
        half->dQuat[0] = gama(half->quat[0], half->AngSpd[0]);
    }
}

}  // namespace

void Simulation::IntegrateDem() {
    const double halfDt = 0.5 * ode.StepSize;
    CudaDemSolver solver;
    solver.initialize(pt, nullptr, control, mechPP, mechPT, mechPW, ode);
    restoreContactHistoryIfPresent(solver);
    const AsyncOutputOptions options = outputOptions(control);
    std::unique_ptr<AsyncOutputWriter> output;
    if (options.any(0)) output = std::make_unique<AsyncOutputWriter>(pt, 0, options);

    std::cout << "CUDA DEM solver: " << CudaDeviceDescription() << std::endl;
    solver.computeForces(0.0, nullptr, nullptr, true);
    std::cout << "Persistent CUDA memory = "
              << solver.deviceMemoryBytes() / (1024.0 * 1024.0) << " MiB" << std::endl;
    enqueueCudaFrame(output.get(), solver, nullptr, nullptr, nullptr, -1, 0.0);
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
            enqueueCudaFrame(output.get(), solver, nullptr, nullptr, nullptr,
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
    saveContactHistory(solver);
}

void Simulation::IntegrateDemMultiBody() {
    auto prescribedMotion = loadPrescribedMotion(
        control.PrescribedMotionFile, ode.StartStep, ode.EndStep);
    std::size_t prescribedMotionCursor = 0;
    applyPrescribedMotion(prescribedMotion, prescribedMotionCursor,
                          ode.StartStep, 0.5 * ode.StepSize, bodyset, nullptr);
    bodyset.updateMatrix();
    const double halfDt = 0.5 * ode.StepSize;
    BodyForceStorage bodyForceStorage(bodyset.Num);
    BodyForceStorage bodyImpulseStorage(bodyset.Num);
    BodyAccelerationStorage bodyAccelerationStorage(bodyset.Num);
    double** bodyForces = bodyForceStorage.data();
    double** bodyImpulses = bodyImpulseStorage.data();
    double** bodyAccelerations = bodyAccelerationStorage.data();
    VHALF bodyHalf(bodyset.Num);
    ACC bodyAcceleration(bodyset.Num);

    CudaDemSolver solver;
    solver.initialize(pt, &bodyset, control, mechPP, mechPT, mechPW, ode);
    restoreContactHistoryIfPresent(solver);
    const AsyncOutputOptions options = outputOptions(control);
    std::unique_ptr<AsyncOutputWriter> output;
    if (options.any(bodyset.Num))
        output = std::make_unique<AsyncOutputWriter>(pt, bodyset.Num, options);

    std::cout << "CUDA DEM-MBD solver: " << CudaDeviceDescription() << std::endl;
    solver.computeForces(0.0, &bodyset, bodyForces, true);
    std::cout << "Persistent CUDA memory = "
              << solver.deviceMemoryBytes() / (1024.0 * 1024.0) << " MiB" << std::endl;
    bodyImpulseStorage.clear();
    enqueueCudaFrame(
        output.get(), solver, &bodyset, bodyForces, bodyImpulses, -1, 0.0);
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
        applyPrescribedMotion(prescribedMotion, prescribedMotionCursor,
                              step, halfDt, bodyset, &bodyHalf);
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
        bodyImpulseStorage.add(bodyForces, ode.StepSize);
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
                             bodyImpulses, step, currentTime);
            bodyImpulseStorage.clear();
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
    saveContactHistory(solver);
}
