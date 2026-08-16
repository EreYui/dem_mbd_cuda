#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

std::vector<double> row(const fs::path& file, int requestedRow = 0) {
    std::ifstream input(file);
    if (!input) throw std::runtime_error("cannot open " + file.string());
    std::string line;
    for (int i = 0; i <= requestedRow; ++i)
        if (!std::getline(input, line)) throw std::runtime_error("missing row in " + file.string());
    std::istringstream values(line);
    std::vector<double> result; double value;
    while (values >> value) result.push_back(value);
    return result;
}

void near(double actual, double expected, double tolerance, const std::string& name) {
    const double error = std::abs(actual - expected);
    std::cout << name << ": actual=" << actual << " expected=" << expected
              << " abs_error=" << error << " tolerance=" << tolerance << '\n';
    if (error > tolerance) throw std::runtime_error(name + " exceeds tolerance");
}

int main(int argc, char** argv) {
    if (argc != 3) return EXIT_FAILURE;
    try {
        const std::string scenario = argv[1];
        const fs::path out = argv[2];
        if (scenario == "free_fall") {
            const auto state = row(out / "state_particles" / "ph.00000010.bt");
            const auto force = row(out / "force_particles" / "F.00000000.bt");
            near(state.at(6), 1.0 + 0.5 * -9.8 * 0.01 * 0.01, 2e-12, "free-fall z");
            near(state.at(9), -9.8 * 0.01, 2e-12, "free-fall vz");
            near(force.at(14), -9.8, 2e-12, "free-fall gravity force");
        } else if (scenario == "wall_contact") {
            const auto force = row(out / "force_particles" / "F.00000000.bt");
            near(force.at(8), 100.0, 1e-9, "wall normal force");
            near(force.at(17), 100.0, 1e-9, "wall resultant force");
        } else if (scenario == "two_particle") {
            const auto first = row(out / "force_particles" / "F.00000000.bt", 0);
            const auto second = row(out / "force_particles" / "F.00000000.bt", 1);
            near(first.at(0), -100.0, 1e-9, "particle 10 contact force x");
            near(second.at(0), 100.0, 1e-9, "particle 20 contact force x");
            near(first.at(0) + second.at(0), 0.0, 1e-12, "particle action-reaction");
        } else if (scenario == "triangle_contact") {
            const auto particle = row(out / "force_particles" / "F.00000000.bt");
            const auto body = row(out / "state_bodys" / "At.00000000.bt");
            near(particle.at(23), 100.0, 1e-8, "triangle force on particle z");
            near(body.at(16), -100.0, 1e-8, "triangle force on body z");
            near(particle.at(23) + body.at(16), 0.0, 1e-10, "triangle action-reaction");
        } else if (scenario == "body_states") {
            const auto fixed = row(out / "state_bodys" / "At.00000010.bt", 0);
            const auto prescribed = row(out / "state_bodys" / "At.00000010.bt", 1);
            near(fixed.at(1), 0.0, 1e-12, "fixed body x");
            near(fixed.at(4), 0.0, 1e-12, "fixed body vx");
            near(prescribed.at(1), 10.02, 1e-11, "prescribed body x");
            near(prescribed.at(4), 2.0, 1e-12, "prescribed body vx");
            near(prescribed.at(9), 1.0, 1e-12, "prescribed body omega z");
            near(prescribed.at(13), std::sin(0.005), 2e-9,
                 "prescribed body quaternion z");
        } else {
            throw std::runtime_error("unknown scenario " + scenario);
        }
    } catch (const std::exception& error) {
        std::cerr << "scenario validation failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
