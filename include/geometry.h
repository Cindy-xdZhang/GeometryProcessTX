#ifndef UTILS_GEOMETRY_H
#define UTILS_GEOMETRY_H

#include <cassert>
#include <vector>
#include <algorithm>
#include "utils/camera.h"
#include <Eigen/Eigen>

extern "C"
{
#include "raytri.c"
}


////////////////////////////////////////////////////////////////////////////////
/// geometry utils
////////////////////////////////////////////////////////////////////////////////

namespace geometry
{
    ////////////////////////////////////////////////////////////////////////////
    /// Intersections
    ////////////////////////////////////////////////////////////////////////////

    struct Hit
    {
        int id = -1; // primitive id
        int gid = -1; // geometry id
        float u = 0.0; // barycentric coordinates
        float v = 0.0;
        float t = 0.0; // distance = direction*t to intersection
    };

    template<
            typename DerivedPos,
            typename DerivedModel,
            typename DerivedProj,
            typename DerivedViewport,
            typename DerivedS,
            typename DerivedDir>
    inline void unproject_ray(
            const Eigen::MatrixBase<DerivedPos>& pos,
            const Eigen::MatrixBase<DerivedModel>& model,
            const Eigen::MatrixBase<DerivedProj>& proj,
            const Eigen::MatrixBase<DerivedViewport>& viewport,
            Eigen::PlainObjectBase<DerivedS>& s,
            Eigen::PlainObjectBase<DerivedDir>& dir)
    {
        // Source and direction on screen
        typedef Eigen::Matrix<typename DerivedS::Scalar, 3, 1> Vec3;
        Vec3 win_s(pos(0, 0), pos(1, 0), 0);
        Vec3 win_d(pos(0, 0), pos(1, 0), 1);
        // Source, destination and direction in world
        Vec3 d;
        utils::unproject(win_s, model, proj, viewport, s);
        utils::unproject(win_d, model, proj, viewport, d);
        dir = d - s;
    }

    template<
            typename DerivedSource,
            typename DerivedDir,
            typename DerivedV,
            typename DerivedF>
    inline bool ray_mesh_intersect(
            const Eigen::MatrixBase<DerivedSource>& s,
            const Eigen::MatrixBase<DerivedDir>& dir,
            const Eigen::MatrixBase<DerivedV>& V,
            const Eigen::MatrixBase<DerivedF>& F,
            std::vector<Hit>& hits)
    {
        // Should be but can't be const
        Eigen::Vector3d s_d = s.template cast<double>();
        Eigen::Vector3d dir_d = dir.template cast<double>();
        hits.clear();
        hits.reserve(F.rows());

        // loop over all triangles
        for (int f = 0; f < F.rows(); f++)
        {
            // Should be but can't be const
            Eigen::RowVector3d v0 = V.row(F(f, 0)).template cast<double>();
            Eigen::RowVector3d v1 = V.row(F(f, 1)).template cast<double>();
            Eigen::RowVector3d v2 = V.row(F(f, 2)).template cast<double>();
            // shoot ray, record hit
            double t, u, v;
            if (intersect_triangle1(
                    s_d.data(), dir_d.data(), v0.data(), v1.data(), v2.data(), &t, &u, &v) &&
                t > 0)
            {
                hits.push_back({(int) f, (int) -1, (float) u, (float) v, (float) t});
            }
        }
        // Sort hits based on distance
        std::sort(
                hits.begin(),
                hits.end(),
                [](const Hit& a, const Hit& b) -> bool
                {
                    return a.t < b.t;
                });
        return !hits.empty();
    }

    template<
            typename DerivedSource,
            typename DerivedDir,
            typename DerivedV,
            typename DerivedF>
    inline bool ray_mesh_intersect(
            const Eigen::MatrixBase<DerivedSource>& source,
            const Eigen::MatrixBase<DerivedDir>& dir,
            const Eigen::MatrixBase<DerivedV>& V,
            const Eigen::MatrixBase<DerivedF>& F,
            Hit& hit)
    {
        std::vector<Hit> hits;
        ray_mesh_intersect(source, dir, V, F, hits);
        if (!hits.empty())
        {
            hit = hits.front();
            return true;
        }
        else
        {
            return false;
        }
    }

    template<typename DerivedObj>
    inline int unproject_in_mesh(
            const Eigen::Vector2f& pos,
            const Eigen::Matrix4f& model,
            const Eigen::Matrix4f& proj,
            const Eigen::Vector4f& viewport,
            const std::function<
                    void(
                            const Eigen::Vector3f&,
                            const Eigen::Vector3f&,
                            std::vector<Hit>&)
            >& shoot_ray,
            Eigen::PlainObjectBase<DerivedObj>& obj,
            std::vector<Hit>& hits)
    {
        Eigen::Vector3f s, dir;
        unproject_ray(pos, model, proj, viewport, s, dir);
        shoot_ray(s, dir, hits);
        switch (hits.size())
        {
            case 0:
                break;
            case 1:
            {
                obj = (s + dir * hits[0].t).cast<typename DerivedObj::Scalar>();
                break;
            }
            case 2:
            default:
            {
                obj = 0.5 * ((s + dir * hits[0].t) + (s + dir * hits[1].t)).cast<typename DerivedObj::Scalar>();
                break;
            }
        }
        return (int) hits.size();
    }

    template<typename DerivedV, typename DerivedF, typename DerivedObj>
    inline int unproject_in_mesh(
            const Eigen::Vector2f& pos,
            const Eigen::Matrix4f& model,
            const Eigen::Matrix4f& proj,
            const Eigen::Vector4f& viewport,
            const Eigen::MatrixBase<DerivedV>& V,
            const Eigen::MatrixBase<DerivedF>& F,
            Eigen::PlainObjectBase<DerivedObj>& obj,
            std::vector<Hit>& hits)
    {
        const auto& shoot_ray = [&V, &F](
                const Eigen::Vector3f& s,
                const Eigen::Vector3f& dir,
                std::vector<Hit>& hits)
        {
            ray_mesh_intersect(s, dir, V, F, hits);
        };
        return unproject_in_mesh(pos, model, proj, viewport, shoot_ray, obj, hits);
    }

    template<typename DerivedV, typename DerivedF, typename DerivedObj>
    inline int unproject_in_mesh(
            const Eigen::Vector2f& pos,
            const Eigen::Matrix4f& model,
            const Eigen::Matrix4f& proj,
            const Eigen::Vector4f& viewport,
            const Eigen::MatrixBase<DerivedV>& V,
            const Eigen::MatrixBase<DerivedF>& F,
            Eigen::PlainObjectBase<DerivedObj>& obj)
    {
        std::vector<Hit> hits;
        return unproject_in_mesh(pos, model, proj, viewport, V, F, obj, hits);
    }
}

#endif // UTILS_GEOMETRY_H
