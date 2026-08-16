#pragma once

#include "matrixtrans.h"

namespace contact_geometry {

// Returns the area of the portion of triangle (a, b, c) that lies inside
// the sphere. Geometrically this is the intersection between the triangle
// and the disk produced by cutting the sphere with the triangle plane.
double sphereTriangleIntersectionArea(
    const vector3d& sphereCenter,
    double sphereRadius,
    const vector3d& a,
    const vector3d& b,
    const vector3d& c);

}  // namespace contact_geometry
