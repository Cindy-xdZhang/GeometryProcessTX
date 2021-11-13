#ifndef VIEWER_H
#define VIEWER_H

#include "viewer/opengl/MeshGL.h"
#include "ViewerCore.h"
#include "ViewerData.h"
#include "ViewerPlugin.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <vector>
#include <string>
#include <cstdint>


struct GLFWwindow;
struct GLFWmonitor;

class Viewer
{

public:

    // UI Enumerations
    enum class MouseButton
    {
        Left, Middle, Right
    };
    enum class MouseMode
    {
        None, Rotation, Zoom, Pan, Translation
    };

    [[maybe_unused]] int launch(
            bool resizable = true, bool fullscreen = false, bool maximize = false,
            const std::string& name = "viewer", int width = 0, int height = 0);

    [[maybe_unused]] int launch_init(
            bool resizable = true, bool fullscreen = false, bool maximize = false,
            const std::string& name = "viewer", int width = 0, int height = 0);

    void launch_rendering(bool loop = true);
    void launch_shut();

    void init();

    void init_plugins();
    void shutdown_plugins();

    Viewer();
    ~Viewer();

    // Mesh IO
    bool load(const std::string& filename, bool camera = true);
    bool save(const std::string& filename);
    bool unload();

    // Screenshots
    int screenshot(const char* filename) const;
    [[nodiscard]] int screenshot(const std::string& filename) const;

    // Callbacks
    bool key_pressed(unsigned int key, int modifier);
    bool key_down(int key, int modifier);
    bool key_up(int key, int modifier);
    bool key_repeat(int key, int modifier);
    bool mouse_down(MouseButton button, int modifier);
    bool mouse_up(MouseButton button, int modifier);
    bool mouse_move(int mouse_x, int mouse_y);
    bool mouse_scroll(float delta_y);

    // Draw everything
    void draw(bool first);
    void refresh();

    // OpenGL context resize
    void resize(int w, int h); // explicitly set window size
    void post_resize(int w, int h); // external resize due to user interaction

    // Helper functions
    void snap_to_canonical_quaternion();
    bool open_dialog_load_mesh();
    bool open_dialog_save_mesh();

    ////////////////////////////////////////////////////////////////////////////
    /// Multi-mesh methods
    ////////////////////////////////////////////////////////////////////////////

    // Return the current mesh, or the mesh corresponding to a given unique identifier
    //
    // Inputs:
    //   mesh_id : unique identifier associated to the desired mesh (current
    //             mesh if -1)
    ViewerData& data(int data_id = -1);
    [[nodiscard]] const ViewerData& data(int data_id = -1) const;

    // Append a new "slot" for a mesh (i.e., create empty entries at the end of
    // the data_list and opengl_state_list).
    //
    // Inputs:
    //   visible : If true, the new mesh is set to be visible on all existing
    //             viewports
    // Returns the id of the last appended mesh
    //
    // Side Effects:
    //   selected_data_index is set this newly created, last entry (i.e.,
    //   #meshes-1)
    int append_data(bool visible = true);

    // Erase a mesh (i.e., its corresponding data and state entries in data_list
    // and opengl_state_list)
    //
    // Inputs:
    //   index : index of mesh to erase
    // Returns whether erasure was successful <=> cannot erase last mesh
    //
    // Side Effects:
    //   If selected_data_index is greater than or equal to index then it is
    //   decremented
    // Example:
    //   // Erase all mesh slots except first and clear remaining mesh
    //   viewer.selected_data_index = viewer.data_list.size()-1;
    //   while(viewer.erase_data(viewer.selected_data_index)){};
    //   viewer.data().clear();
    //
    bool erase_data(size_t index);

    ////////////////////////////////////////////////////////////////////////////
    /// Multi-viewport methods
    ////////////////////////////////////////////////////////////////////////////

    // Return the current viewport, or the viewport corresponding to a given unique identifier
    //
    // Inputs:
    //   core_id : unique identifier corresponding to the desired viewport
    //             (current viewport if 0)
    ViewerCore& core(unsigned core_id = 0);
    [[nodiscard]] const ViewerCore& core(unsigned core_id = 0) const;

    // Append a new "slot" for a viewport (i.e., copy properties of the current viewport, only
    // changing the viewport size/position)
    //
    // Inputs:
    //   viewport     : Vector specifying the viewport origin and size in screen
    //                  coordinates.
    //   append_empty : If true, existing meshes are hidden on the new viewport.
    //
    // Returns the unique id of the newly inserted viewport. There can be a maximum of 31
    //   viewports created in the same viewport. Erasing a viewport does not change the id of
    //   other existing viewports
    int append_core(Eigen::Vector4f viewport, bool append_empty = false);

    // Erase a viewport
    //
    // Inputs:
    //   index : index of the viewport to erase
    bool erase_core(size_t index);

    // Retrieve viewport index from its unique identifier
    // Returns 0 if not found
    [[nodiscard]] size_t core_index(int id) const;

    // Change selected_core_index to the viewport containing the mouse
    // (current_mouse_x, current_mouse_y)
    void select_hovered_core();

public:

    ////////////////////////////////////////////////////////////////////////////
    /// Member variables
    ////////////////////////////////////////////////////////////////////////////

    std::string path = ".";

    GLFWwindow* window;

    // Stores all the meshes that should be visualized
    std::vector<ViewerData> data_list;
    int selected_data_index;
    int next_data_id;

    // Stores all the viewing options
    std::vector<ViewerCore> core_list;
    int selected_core_index;
    int next_core_id;

    // List of registered plugins
    std::vector<ViewerPlugin*> plugins;

    // Temporary data stored when the mouse button is pressed
    MouseMode mouse_mode = Viewer::MouseMode::None;
    Eigen::Quaternionf down_rotation;
    int current_mouse_x = 0;
    int current_mouse_y = 0;
    int down_mouse_x = 0;
    int down_mouse_y = 0;
    float down_mouse_z = 0;
    Eigen::Vector3f down_translation;
    bool down;
    bool hack_never_moved;

    // Keep track of the global position of the scroll wheel
    float scroll_position;

    // C++-style functions
    //
    // Returns **true** if action should be cancelled.
    std::function<bool(Viewer& viewer)> callback_init;
    std::function<bool(Viewer& viewer)> callback_pre_draw;
    std::function<bool(Viewer& viewer)> callback_post_draw;
    std::function<bool(Viewer& viewer, int w, int h)> callback_post_resize;
    std::function<bool(Viewer& viewer, int button, int modifier)> callback_mouse_down;
    std::function<bool(Viewer& viewer, int button, int modifier)> callback_mouse_up;
    std::function<bool(Viewer& viewer, int mouse_x, int mouse_y)> callback_mouse_move;
    std::function<bool(Viewer& viewer, float delta_y)> callback_mouse_scroll;
    std::function<bool(Viewer& viewer, unsigned int key, int modifiers)> callback_key_pressed;
    std::function<bool(Viewer& viewer, unsigned int key, int modifiers)> callback_key_down;
    std::function<bool(Viewer& viewer, unsigned int key, int modifiers)> callback_key_up;
    std::function<bool(Viewer& viewer, unsigned int key, int modifiers)> callback_key_repeat;

    // Pointers to per-callback data
    void* callback_init_data;
    void* callback_pre_draw_data;
    void* callback_post_draw_data;
    void* callback_mouse_down_data;
    void* callback_mouse_up_data;
    void* callback_mouse_move_data;
    void* callback_mouse_scroll_data;
    void* callback_key_pressed_data;
    void* callback_key_down_data;
    void* callback_key_up_data;
    void* callback_key_repeat_data;

public:
    //use Eigen internal operator new,delete,new[],delete[]
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif // VIEWER_H
