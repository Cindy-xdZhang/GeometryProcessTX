#ifndef UTILS_CONSTANTS_H
#define UTILS_CONSTANTS_H

#include "utils/OpenMesh.h"


namespace utils
{

    // pi
#ifdef M_PI
    const double PI = M_PI;
    const double PI_2 = M_PI_2;
#else
    const double PI = 3.1415926535897932384626433832795;
    const double PI_2 = 1.57079632679489661923;
#endif
    const double PI_3_2 = 4.71238898038468985769; // 3*pi/2

    // Define a standard value for double epsilon
    const double DOUBLE_EPS = 1.0e-14;
    const double DOUBLE_EPS_SQ = 1.0e-28;
    const float FLOAT_EPS = 1.0e-7f;
    const float FLOAT_EPS_SQ = 1.0e-14f;

    // Function returning EPS for corresponding type
    template<typename S_type>
    inline S_type EPS();
    template<typename S_type>
    inline S_type EPS_SQ();

    // Template specializations for float and double
    template<>
    inline float EPS()
    {
        return FLOAT_EPS;
    }

    template<>
    inline double EPS()
    {
        return DOUBLE_EPS;
    }

    template<>
    inline float EPS_SQ()
    {
        return FLOAT_EPS_SQ;
    }

    template<>
    inline double EPS_SQ()
    {
        return DOUBLE_EPS_SQ;
    }

    // axis
    const OpenMesh::Mesh::Point X_AXIS = OpenMesh::Mesh::Point(1.0, 0.0, 0.0);
    const OpenMesh::Mesh::Point Y_AXIS = OpenMesh::Mesh::Point(0.0, 1.0, 0.0);
    const OpenMesh::Mesh::Point Z_AXIS = OpenMesh::Mesh::Point(0.0, 0.0, 1.0);
}

#endif // UTILS_CONSTANTS_H
