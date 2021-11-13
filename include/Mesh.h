#ifndef MESH_H
#define MESH_H

#include <Eigen/Eigen>
#include "utils/colormap.h"
#include "utils/OpenMesh.h"


class ViewerData;

class Mesh
{
    static size_t counter;

public: // enums

    enum NoiseDirection
    {
        NORMAL,
        RANDOM
    };

public:

    Mesh();

    bool load(const std::string& filename);
    bool save(const std::string& filename);

    [[nodiscard]] std::string name() const;
    [[nodiscard]] unsigned int id() const;

    void updateViewerData(ViewerData& data) const;

    OpenMesh::Mesh& mesh();
    [[nodiscard]] OpenMesh::Mesh mesh() const;

    void clean();
    void resetMesh();
    void setMesh(const OpenMesh::Mesh& mesh);

    bool& renderFlatFaces();
    bool& renderEdges();
    Eigen::Vector4f& edgeColor();
    colormap::ColorMapType& colormap();

    [[nodiscard]] size_t numVertices() const;
    [[nodiscard]] Eigen::Vector3d vertex(unsigned int index) const;
    void setVertex(unsigned int index, const Eigen::Vector3d& v);

    [[nodiscard]] Eigen::Vector3d faceCenter(unsigned int index) const;
    void setFaceCenter(unsigned int index, const Eigen::Vector3d& v);

    void normalize();
    double averageEdgeLength();
    double averageDihedralAngle();
    void noise(double standardDeviation = 0.1, NoiseDirection noiseDirection = NORMAL);

    void LaplacianSmoothing(double lambda, int iterations, bool uniform);
    void OptimizingSmoothing(double lambda, double mu, double gama, double theta, int iterations = 5, bool uniform = false);
private:

    OpenMesh::Mesh mMesh;
    OpenMesh::Mesh mMeshOriginal;

    // mesh properties
    unsigned int mID;
    std::string mName;

    // rendering settings
    bool mRenderFlatFaces;
    bool mRenderEdges;
    Eigen::Vector4f mEdgeColor;
    colormap::ColorMapType mColorMap;
};


#endif // MESH_H
