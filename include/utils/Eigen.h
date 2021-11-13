#ifndef UTILS_EIGEN_H
#define UTILS_EIGEN_H

#include <random>
#include <iostream>
#include <Eigen/Eigen>


////////////////////////////////////////////////////////////////////////////////
/// Eigen utils
////////////////////////////////////////////////////////////////////////////////

namespace Eigen
{
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrixXd;

    inline Eigen::Vector3d randomVector()
    {
        std::random_device device;
        std::mt19937 mt(device());
        std::uniform_real_distribution<double> dist(-1.0, 1.0);

        double x, y, z, length;
        do
        {
            x = dist(mt);
            y = dist(mt);
            length = x * x + y * y;
        }
        while (length > 1);

        const double r = 2 * std::sqrt(1 - length);

        x *= r;
        y *= r;
        z = 1 - 2 * length;

        return {x, y, z};
    }
}

#endif // UTILS_EIGEN_H
