#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct Metrics {
    std::size_t values = 0;
    double maxAbs = 0.0;
    double maxRel = 0.0;
    std::size_t maxAbsIndex = 0;
    std::size_t nonFinite = 0;
};

std::vector<double> readNumbers(const fs::path& file) {
    std::ifstream input(file);
    if (!input) throw std::runtime_error("cannot open " + file.string());
    std::vector<double> values;
    double value = 0.0;
    while (input >> value) values.push_back(value);
    if (!input.eof()) throw std::runtime_error("non-numeric token in " + file.string());
    return values;
}

Metrics compare(const std::vector<double>& reference,
                const std::vector<double>& candidate) {
    if (reference.size() != candidate.size()) {
        throw std::runtime_error("value count differs: " +
            std::to_string(reference.size()) + " vs " +
            std::to_string(candidate.size()));
    }
    Metrics result;
    result.values = reference.size();
    for (std::size_t i = 0; i < reference.size(); ++i) {
        if (!std::isfinite(reference[i]) || !std::isfinite(candidate[i])) {
            ++result.nonFinite;
            continue;
        }
        const double absolute = std::abs(candidate[i] - reference[i]);
        const double scale = std::max({std::abs(reference[i]),
                                       std::abs(candidate[i]), 1.0e-12});
        const double relative = absolute / scale;
        if (absolute > result.maxAbs) {
            result.maxAbs = absolute;
            result.maxAbsIndex = i;
        }
        result.maxRel = std::max(result.maxRel, relative);
    }
    return result;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: output_compare REFERENCE_OUTPUT CANDIDATE_OUTPUT\n";
        return EXIT_FAILURE;
    }
    try {
        const fs::path referenceRoot = fs::absolute(argv[1]);
        const fs::path candidateRoot = fs::absolute(argv[2]);
        std::cout << "file,values,max_abs,max_rel,max_abs_index,non_finite\n";
        for (const auto& entry : fs::recursive_directory_iterator(referenceRoot)) {
            if (!entry.is_regular_file() || entry.path().extension() != ".bt") continue;
            const fs::path relative = fs::relative(entry.path(), referenceRoot);
            const fs::path candidate = candidateRoot / relative;
            if (!fs::is_regular_file(candidate)) {
                throw std::runtime_error("missing candidate file " + candidate.string());
            }
            const auto referenceValues = readNumbers(entry.path());
            const auto candidateValues = readNumbers(candidate);
            const Metrics metrics = compare(referenceValues, candidateValues);
            std::cout << relative.generic_string() << ',' << metrics.values << ','
                      << std::scientific << std::setprecision(8)
                      << metrics.maxAbs << ',' << metrics.maxRel << ','
                      << metrics.maxAbsIndex << ',' << metrics.nonFinite << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "output comparison failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
