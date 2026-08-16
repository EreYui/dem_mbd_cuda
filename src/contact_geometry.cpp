#include "contact_geometry.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace contact_geometry {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kTolerance = 64.0 * std::numeric_limits<double>::epsilon();

double signedCross(
    const vector3d& a,
    const vector3d& b,
    const vector3d& unitNormal) {
    return unitNormal.dot(a.cross(b));
}

double sectorArea(
    const vector3d& a,
    const vector3d& b,
    const vector3d& unitNormal) {
    return 0.5 * std::atan2(signedCross(a, b, unitNormal), a.dot(b));
}

// Signed area contribution of one polygon edge to the intersection between
// the polygon and a unit disk centered at the origin. The edge is split at
// its circle intersections. Chords inside the disk contribute triangle area;
// portions outside the disk contribute circular-sector area.
double unitDiskEdgeArea(
    const vector3d& a,
    const vector3d& b,
    const vector3d& unitNormal) {
    const vector3d direction = b - a;
    const double quadraticA = direction.squaredNorm();
    if (quadraticA <= kTolerance) {
        return 0.0;
    }

    const double quadraticB = 2.0 * a.dot(direction);
    const double quadraticC = a.squaredNorm() - 1.0;
    double discriminant = quadraticB * quadraticB
        - 4.0 * quadraticA * quadraticC;
    const double discriminantTolerance = kTolerance
        * (quadraticB * quadraticB
            + std::abs(4.0 * quadraticA * quadraticC) + 1.0);

    std::array<double, 4> parameters{{0.0, 1.0, 0.0, 0.0}};
    std::size_t parameterCount = 2;

    if (discriminant >= -discriminantTolerance) {
        discriminant = std::max(0.0, discriminant);
        const double squareRoot = std::sqrt(discriminant);

        // The parameters are only used to split a finite segment and are
        // clamped by the interval checks below.
        const double denominator = 2.0 * quadraticA;
        const double roots[2] = {
            (-quadraticB - squareRoot) / denominator,
            (-quadraticB + squareRoot) / denominator
        };

        for (double root : roots) {
            if (root > kTolerance && root < 1.0 - kTolerance) {
                parameters[parameterCount++] = root;
            }
        }
    }

    std::sort(parameters.begin(), parameters.begin() + parameterCount);
    const auto uniqueEnd = std::unique(
        parameters.begin(), parameters.begin() + parameterCount,
        [](double lhs, double rhs) {
            return std::abs(lhs - rhs) <= kTolerance;
        });
    parameterCount = static_cast<std::size_t>(
        std::distance(parameters.begin(), uniqueEnd));

    double area = 0.0;
    for (std::size_t i = 0; i + 1 < parameterCount; ++i) {
        const double t0 = parameters[i];
        const double t1 = parameters[i + 1];
        const vector3d p = a + direction * t0;
        const vector3d q = a + direction * t1;
        const vector3d midpoint = a + direction * (0.5 * (t0 + t1));

        if (midpoint.squaredNorm() <= 1.0 + kTolerance) {
            area += 0.5 * signedCross(p, q, unitNormal);
        } else {
            area += sectorArea(p, q, unitNormal);
        }
    }

    return area;
}

}  // namespace

double sphereTriangleIntersectionArea(
    const vector3d& sphereCenter,
    double sphereRadius,
    const vector3d& a,
    const vector3d& b,
    const vector3d& c) {
    if (!(sphereRadius > 0.0) || !std::isfinite(sphereRadius)) {
        return 0.0;
    }

    const vector3d unnormalizedNormal = (b - a).cross(c - a);
    const double doubleTriangleArea = unnormalizedNormal.norm();
    if (!(doubleTriangleArea > 0.0) || !std::isfinite(doubleTriangleArea)) {
        return 0.0;
    }

    const vector3d unitNormal = unnormalizedNormal / doubleTriangleArea;
    const double signedDistance = unitNormal.dot(sphereCenter - a);
    const double sphereRadiusSquared = sphereRadius * sphereRadius;
    const double diskRadiusSquared = sphereRadiusSquared
        - signedDistance * signedDistance;

    if (!(diskRadiusSquared > 0.0) || !std::isfinite(diskRadiusSquared)) {
        return 0.0;
    }

    const double diskRadius = std::sqrt(diskRadiusSquared);
    const vector3d diskCenter = sphereCenter - unitNormal * signedDistance;

    // Normalize the disk to unit radius. This keeps tolerances independent of
    // the physical scale and mirrors the conditioning strategy of the former
    // general-purpose overlap implementation.
    const vector3d normalizedA = (a - diskCenter) / diskRadius;
    const vector3d normalizedB = (b - diskCenter) / diskRadius;
    const vector3d normalizedC = (c - diskCenter) / diskRadius;

    double normalizedArea = 0.0;
    normalizedArea += unitDiskEdgeArea(normalizedA, normalizedB, unitNormal);
    normalizedArea += unitDiskEdgeArea(normalizedB, normalizedC, unitNormal);
    normalizedArea += unitDiskEdgeArea(normalizedC, normalizedA, unitNormal);
    normalizedArea = std::abs(normalizedArea);

    const double triangleArea = 0.5 * doubleTriangleArea;
    const double diskArea = kPi * diskRadiusSquared;
    const double maximumArea = std::min(triangleArea, diskArea);
    const double area = normalizedArea * diskRadiusSquared;

    if (!std::isfinite(area)) {
        return 0.0;
    }
    return std::clamp(area, 0.0, maximumArea);
}

}  // namespace contact_geometry
