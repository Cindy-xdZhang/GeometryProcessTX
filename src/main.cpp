#include <filesystem>

#include <Eigen/src/Core/util/DisableStupidWarnings.h>
#include <Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Vector2d)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Vector3d)

#include "viewer/Viewer.h"
#include "viewer/plugins/WidgetsPlugin.h"
#include "viewer/plugins/LoaderPlugin.h"

#include "viewer/widgets/MainWidget.h"
#include "viewer/widgets/MeshWidget.h"


int main(int /*argc*/, char* /*argv*/[])
{
    Viewer viewer;

    // Change default path
    viewer.path = "/path/to/workspace/";

    // Attach menu plugins
    LoaderPlugin loader;
    viewer.plugins.push_back(&loader);
    WidgetsPlugin menus;
    viewer.plugins.push_back(&menus);
    menus.add(new MainWidget(true, &loader));
    menus.add(new MeshWidget(true, &loader));

    // General drawing settings
    viewer.core().is_animating = false;
    viewer.core().animation_max_fps = 30.;
    viewer.core().background_color = Eigen::Vector4f(0.6f, 0.6f, 0.6f, 1.0f);

    // Initialize viewer
    if (viewer.launch_init(true, false, true, "viewer", 0, 0) == EXIT_FAILURE)
    {
        viewer.launch_shut();
        return EXIT_FAILURE;
    }

    // Example of loading a model programmatically
#ifdef PROJECT_DIRECTORY
    std::filesystem::path path =
            std::filesystem::path(PROJECT_DIRECTORY)
                    .append("models").append("botanic_garden.obj");
    viewer.load(path.string());
#endif

    // Rendering
    viewer.launch_rendering(true);
    viewer.launch_shut();
    return EXIT_SUCCESS;
}
