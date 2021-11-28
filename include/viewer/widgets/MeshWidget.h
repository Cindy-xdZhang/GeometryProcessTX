#ifndef MESH_WIDGET_H
#define MESH_WIDGET_H
//#include <OpenMesh/Tools/Subdivider/Uniform/CatmullClarkT.hh>
#include "viewer/ViewerWidget.h"
#include "viewer/plugins/LoaderPlugin.h"
#include "geometry.h"

class MeshWidget : public ViewerWidget
{

public:

    explicit MeshWidget(bool expanded, LoaderPlugin* loader = nullptr);

    ImVec4 draw(bool first, float scaling, float xSize, float ySize, float xPos, float yPos) override;
	bool mouse_down(int button, int modifier) override;
	//bool mouse_up(int button, int modifier) override;
    bool mouse_move(int mouse_x, int mouse_y) override;
    
    Eigen::Vector3d pick(int& face_index, int& vertex_index);
    Eigen::Vector3d pick(int& face_index, int& vertex_index, geometry::Hit& hit);
private:

    LoaderPlugin* mLoader;

    double noiseStandardDeviation = 0.1;
    double noiseStandardDeviationMin = 0.0;
    double noiseStandardDeviationMax = 5.0;
    Mesh::NoiseDirection noiseDirection = Mesh::NoiseDirection::NORMAL;
    int current_mouse_x;
    int current_mouse_y;
    int down_mouse_x;
    int down_mouse_y;
    int picked_face_index;
    int picked_vertex_index=-1;
};

#endif // MESH_WIDGET_H
