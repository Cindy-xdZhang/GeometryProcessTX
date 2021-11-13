#ifndef UTILS_MATH_H
#define UTILS_MATH_H

#include <cmath>
#include <random>

#include <Eigen/Eigen>

#include "utils/constants.h"


////////////////////////////////////////////////////////////////////////////////
/// math utils
////////////////////////////////////////////////////////////////////////////////

namespace maths
{
    template<typename T>
    int sgn(T val)
    {
        //return (val < 0.0) ? -1 : 1;
        return std::copysign(1, val);
        //return (T(0) < val) - (val < T(0));
    }

    inline double degree2rad(double degree)
    {
        return degree * utils::PI / 180.0;
    }

    // http://www.antisphere.com/Wiki/tools:anttweakbar

    inline double dot(const double* a, const double* b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    inline void cross(const double* a, const double* b, double* out)
    {
        out[0] = a[1] * b[2] - a[2] * b[1];
        out[1] = a[2] * b[0] - a[0] * b[2];
        out[2] = a[0] * b[1] - a[1] * b[0];
    }

    // https://math.stackexchange.com/a/3582461

    inline Eigen::Vector3d normal(const Eigen::Vector3d& v)
    {
        const double s = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        const double g = std::copysign(s, v[2]);
        const double h = v[2] + g;
        return {g * h - v[0] * v[0], -v[0] * v[1], -v[0] * h};
    }

    // random

    inline double gaussian(double mean, double std)
    {
        std::random_device device;
        std::mt19937 mt(device());
        std::uniform_real_distribution<double> dist(-1.0, 1.0);

        static double v1, v2, s;
        static int phase = 0;
        double x;

        if (phase == 0)
        {
            do
            {
                v1 = dist(mt);
                v2 = dist(mt);
                s = v1 * v1 + v2 * v2;
            }
            while (s >= 1 || s == 0);

            x = v1 * sqrt(-2 * log(s) / s);
        }
        else
        {
            x = v2 * sqrt(-2 * log(s) / s);
        }

        phase = 1 - phase;

        return x * std + mean;
    }
}

#endif // UTILS_MATH_H
