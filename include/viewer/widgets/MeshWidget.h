#ifndef MESH_WIDGET_H
#define MESH_WIDGET_H
//#include <OpenMesh/Tools/Subdivider/Uniform/CatmullClarkT.hh>
#include "viewer/ViewerWidget.h"
#include "viewer/plugins/LoaderPlugin.h"


class MeshWidget : public ViewerWidget
{

public:

    explicit MeshWidget(bool expanded, LoaderPlugin* loader = nullptr);

    ImVec4 draw(bool first, float scaling, float xSize, float ySize, float xPos, float yPos) override;

private:

    LoaderPlugin* mLoader;

    double noiseStandardDeviation = 0.1;
    double noiseStandardDeviationMin = 0.0;
    double noiseStandardDeviationMax = 5.0;
    Mesh::NoiseDirection noiseDirection = Mesh::NoiseDirection::NORMAL;
};

#endif // MESH_WIDGET_H
