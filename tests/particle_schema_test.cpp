#include "particle.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

bool close(double first, double second)
{
    return std::abs(first - second) <= 1.0e-15;
}

void write(const fs::path& path, const std::string& content)
{
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create test fixture " + path.string());
    output << content;
}

bool rejects(PARTICLE& particles, const fs::path& path, bool inertia)
{
    try {
        if (inertia) particles.LoadPrincipalInertia(path.string());
        else particles.LoadClumpComponents(path.string());
    }
    catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

} // namespace

int main()
{
    const fs::path directory = fs::temp_directory_path()
        / ("dem_mbd_particle_schema_" + std::to_string(std::rand()));
    try {
        fs::create_directories(directory);
        const fs::path particleFile = directory / "Particles.bt";
        write(
            particleFile,
            "31 0 2.0 0.1 0 0 0 0 0 0 1 0 0 0 0 0 0\n"
            "47 0 3.0 0.2 1 0 0 0 0 0 1 0 0 0 0 0 0\n");

        PARTICLE particles;
        particles.LoadParticles(particleFile.string());
        particles.LoadPrincipalInertia("");
        particles.LoadClumpComponents("");
        require(close(particles.Inertia[0][0], 0.008), "legacy sphere inertia changed");
        require(particles.ComponentNum == 2, "legacy sphere topology count changed");
        require(particles.ComponentId[1] == 47 && particles.ComponentOwner[1] == 1,
                "legacy sphere topology identity changed");

        const fs::path inertiaFile = directory / "ParticleInertia.bt";
        write(inertiaFile, "31 0.01 0.02 0.03\n47 0.04 0.05 0.06\n");
        particles.LoadPrincipalInertia(inertiaFile.string());
        require(close(particles.Inertia[0][1], 0.02), "explicit inertia was not loaded");
        require(close(particles.Inertia[1][2], 0.06), "explicit inertia owner mapping failed");

        const fs::path componentsFile = directory / "ClumpComponents.bt";
        write(
            componentsFile,
            "101 31 0.08 -0.02 0 0\n"
            "102 31 0.08 0.02 0 0\n"
            "201 47 0.20 0 0 0\n");
        particles.LoadClumpComponents(componentsFile.string());
        require(particles.ComponentNum == 3, "explicit topology count failed");
        require(particles.ComponentOwner[1] == 0 && particles.ComponentId[2] == 201,
                "explicit topology owner mapping failed");

        const fs::path badInertia = directory / "bad_inertia.bt";
        write(badInertia, "31 0.01 0.02 0.03\n47 0.04 -0.05 0.06\n");
        require(rejects(particles, badInertia, true), "negative inertia did not fail closed");

        const fs::path badComponents = directory / "bad_components.bt";
        write(badComponents, "101 31 0.08 0 0 0\n101 47 0.20 0 0 0\n");
        require(rejects(particles, badComponents, false),
                "duplicate component ID did not fail closed");

        const fs::path negativeComponent = directory / "negative_component.bt";
        write(negativeComponent, "-1 31 0.08 0 0 0\n201 47 0.20 0 0 0\n");
        require(rejects(particles, negativeComponent, false),
                "negative component ID did not fail closed");

        const fs::path crowdedComponents = directory / "crowded_components.bt";
        write(
            crowdedComponents,
            "100 31 0.01 0 0 0\n101 31 0.01 0 0 0\n"
            "102 31 0.01 0 0 0\n103 31 0.01 0 0 0\n"
            "104 31 0.01 0 0 0\n105 31 0.01 0 0 0\n"
            "106 31 0.01 0 0 0\n107 31 0.01 0 0 0\n"
            "108 31 0.01 0 0 0\n201 47 0.20 0 0 0\n");
        require(rejects(particles, crowdedComponents, false),
                "component owner capacity did not fail closed");
        fs::remove_all(directory);
    }
    catch (const std::exception& error) {
        fs::remove_all(directory);
        std::cerr << "particle schema test failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    std::cout << "particle schema tests passed\n";
    return EXIT_SUCCESS;
}
