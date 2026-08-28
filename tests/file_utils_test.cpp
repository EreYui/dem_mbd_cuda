#include "file_utils.h"

#include <cstdlib>
#include <iostream>
#include <string>

int main() {
    try {
        file_utils::setCaseDirectory(
            file_utils::normalizePath("validation/cases/free_fall"));
        // Deliberately use legacy Windows separators. This must also resolve
        // on Ubuntu when CTest runs with the project root as working directory.
        const auto existingPath = file_utils::requireInputFile(
            "Data\\DATA\\InputFile\\SettingData\\world.par",
            "cross-platform path test");
        if (existingPath.generic_string()
            != "validation/cases/free_fall/InputFile/SettingData/world.par") {
            std::cerr << "unexpected normalized path: "
                      << existingPath.generic_string() << '\n';
            return EXIT_FAILURE;
        }

        file_utils::setCaseDirectory(
            file_utils::normalizePath("validation/cases/body_states"));
        const auto redirectedLegacyPath = file_utils::resolveCasePath(
            "Data\\DATA\\InputFile\\Particles\\Particles.bt");
        if (redirectedLegacyPath.generic_string()
            != "validation/cases/body_states/InputFile/Particles/Particles.bt") {
            std::cerr << "legacy case path was not redirected: "
                      << redirectedLegacyPath.generic_string() << '\n';
            return EXIT_FAILURE;
        }
        const auto caseRelativePath = file_utils::resolveCasePath(
            "OutputFile/state_particles/ph.00001.bt");
        if (caseRelativePath.generic_string()
            != "validation/cases/body_states/OutputFile/state_particles/ph.00001.bt") {
            std::cerr << "case-relative path was not resolved: "
                      << caseRelativePath.generic_string() << '\n';
            return EXIT_FAILURE;
        }
        file_utils::setCaseDirectory(
            file_utils::normalizePath("validation/cases/free_fall"));

        bool missingFileReported = false;
        try {
            file_utils::requireInputFile(
                "Data\\DATA\\InputFile\\definitely_missing_for_test.bt",
                "missing-file test");
        }
        catch (const std::runtime_error& error) {
            const std::string message = error.what();
            missingFileReported = message.find("File not found")
                != std::string::npos;
        }

        if (!missingFileReported) {
            std::cerr << "missing file did not produce the expected error\n";
            return EXIT_FAILURE;
        }
    }
    catch (const std::exception& error) {
        std::cerr << "file_utils test failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }

    std::cout << "file_utils tests passed\n";
    return EXIT_SUCCESS;
}
