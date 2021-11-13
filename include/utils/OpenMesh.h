#ifndef UTILS_OPEN_MESH_H
#define UTILS_OPEN_MESH_H

#include <Eigen/Eigen>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/EigenVectorT.hh>


////////////////////////////////////////////////////////////////////////////////
/// OpenMesh utils
////////////////////////////////////////////////////////////////////////////////

namespace OpenMesh
{
    struct Traits : public OpenMesh::DefaultTraits
    {
        using Point = Eigen::Vector3d;
        using Normal = Eigen::Vector3d;
        using TexCoord2D = Eigen::Vector2d;

        VertexAttributes(OpenMesh::Attributes::Normal);
        HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
        EdgeAttributes(0);
        FaceAttributes(OpenMesh::Attributes::Normal);
    };

    typedef OpenMesh::PolyMesh_ArrayKernelT<Traits> Mesh;

    inline double dot(const OpenMesh::Mesh::Point& a, const OpenMesh::Mesh::Point& b)
    {
        return a.dot(b);
    }

    inline OpenMesh::Mesh::Point cross(const OpenMesh::Mesh::Point& a, const OpenMesh::Mesh::Point& b)
    {
        return a.cross(b);
    }

    inline double length(const OpenMesh::Mesh::Point& a)
    {
        return a.norm();
    }

    inline double sqrlength(const OpenMesh::Mesh::Point& a)
    {
        return a.squaredNorm();
    }

    inline double norm(const OpenMesh::Mesh::Point& a)
    {
        return a.norm();
    }

    inline double sqrnorm(const OpenMesh::Mesh::Point& a)
    {
        return a.squaredNorm();
    }

    inline OpenMesh::Mesh::Point normalize(const OpenMesh::Mesh::Point& a)
    {
        return a.normalized();
    }

    inline bool coface(const OpenMesh::Mesh::VertexIter& u, const OpenMesh::Mesh::VertexIter& v)
    {
        OpenMesh::SmartFaceHandle f11 = u->halfedge().face();
        OpenMesh::SmartFaceHandle f12 = u->halfedge().opp().face();
        OpenMesh::SmartFaceHandle f21 = v->halfedge().face();
        OpenMesh::SmartFaceHandle f22 = v->halfedge().opp().face();
        return f11 == f21 || f11 == f22 || f12 == f21 || f12 == f22;
    }

    // Rendering

    void toRenderData(
            OpenMesh::Mesh mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F,
            Eigen::MatrixXi& FtoF, Eigen::MatrixXd& N, Eigen::MatrixXd& UV,
            Eigen::MatrixXd& P1, Eigen::MatrixXd& P2);
}

#endif // UTILS_OPEN_MESH_H
