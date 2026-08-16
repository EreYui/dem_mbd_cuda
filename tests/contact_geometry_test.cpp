#include "contact_geometry.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

bool nearlyEqual(double actual, double expected, double relativeTolerance = 1e-10) {
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    return std::abs(actual - expected) <= relativeTolerance * scale;
}

void expectArea(
    const std::string& name,
    double expected,
    const vector3d& center,
    double radius,
    const vector3d& a,
    const vector3d& b,
    const vector3d& c) {
    const double actual = contact_geometry::sphereTriangleIntersectionArea(
        center, radius, a, b, c);
    if (!nearlyEqual(actual, expected)) {
        std::cerr << name << " failed: expected " << expected
                  << ", got " << actual << '\n';
        std::exit(EXIT_FAILURE);
    }
}

}  // namespace

int main() {
    const vector3d origin{0.0, 0.0, 0.0};
    const vector3d largeA{-10.0, -10.0, 0.0};
    const vector3d largeB{10.0, -10.0, 0.0};
    const vector3d largeC{0.0, 10.0, 0.0};

    expectArea("disk inside triangle", kPi, origin, 1.0,
        largeA, largeB, largeC);
    expectArea("offset disk inside triangle", 0.64 * kPi,
        {0.0, 0.0, 0.6}, 1.0, largeA, largeB, largeC);
    expectArea("sphere tangent to plane", 0.0,
        {0.0, 0.0, 1.0}, 1.0, largeA, largeB, largeC);
    expectArea("sphere misses plane", 0.0,
        {0.0, 0.0, 2.0}, 1.0, largeA, largeB, largeC);

    expectArea("triangle inside disk", 0.5, origin, 10.0,
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
    expectArea("disjoint", 0.0, origin, 1.0,
        {2.0, 2.0, 0.0}, {3.0, 2.0, 0.0}, {2.0, 3.0, 0.0});
    expectArea("half disk", 0.5 * kPi, origin, 1.0,
        {0.0, -10.0, 0.0}, {10.0, 0.0, 0.0}, {0.0, 10.0, 0.0});
    expectArea("quarter disk", 0.25 * kPi, origin, 1.0,
        {0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}, {0.0, 10.0, 0.0});
    expectArea("reversed orientation", 0.25 * kPi, origin, 1.0,
        {0.0, 10.0, 0.0}, {10.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    expectArea("degenerate triangle", 0.0, origin, 1.0,
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0});
    expectArea("zero-radius sphere", 0.0, origin, 0.0,
        largeA, largeB, largeC);

    const double scale = 1.0e6;
    expectArea("large-scale invariance", kPi * scale * scale,
        {3.0 * scale, -2.0 * scale, 5.0 * scale}, scale,
        {-7.0 * scale, -12.0 * scale, 5.0 * scale},
        {13.0 * scale, -12.0 * scale, 5.0 * scale},
        {3.0 * scale, 8.0 * scale, 5.0 * scale});

    std::cout << "contact_geometry tests passed\n";
    return EXIT_SUCCESS;
}
