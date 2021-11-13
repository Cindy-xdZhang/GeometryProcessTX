#ifndef VIEWER_PLUGIN_H
#define VIEWER_PLUGIN_H

#include <string>
#include <vector>


class Viewer;

class ViewerPlugin
{

public:

    ViewerPlugin();

    virtual ~ViewerPlugin();

    // This function is called when the viewer is initialized (no mesh will be loaded at this stage)
    virtual void init(Viewer* _viewer);

    // This function is called before shutdown
    virtual void shutdown();

    // This function is called before a mesh is loaded
    virtual bool load(const std::string& filename);

    // This function is called before a mesh is unloaded
    virtual bool unload();

    // This function is called before a mesh is saved
    virtual bool save(const std::string& filename);

    // This function is called when the scene is serialized
    virtual bool serialize(std::vector<char>& buffer) const;

    // This function is called when the scene is deserialized
    virtual bool deserialize(const std::vector<char>& buffer);

    // Runs immediately after a new mesh has been loaded.
    virtual bool post_load();

    // This function is called before the draw procedure of Preview3D
    virtual bool pre_draw(bool first);

    // This function is called after the draw procedure of Preview3D
    virtual bool post_draw(bool first);

    // This function is called after the window has been resized
    virtual bool post_resize(int w, int h);

    // This function is called when the mouse button is pressed
    // - button can be GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON or GLUT_RIGHT_BUTTON
    // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
    virtual bool mouse_down(int button, int modifier);

    // This function is called when the mouse button is released
    // - button can be GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON or GLUT_RIGHT_BUTTON
    // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
    virtual bool mouse_up(int button, int modifier);

    // This function is called every time the mouse is moved
    // - mouse_x and mouse_y are the new coordinates of the mouse pointer in screen coordinates
    virtual bool mouse_move(int mouse_x, int mouse_y);

    // This function is called every time the scroll wheel is moved
    // Note: this callback is not working with every glut implementation
    virtual bool mouse_scroll(float delta_y);

    // This function is called when a keyboard key is pressed. Unlike key_down
    // this will reveal the actual character being sent (not just the physical
    // key)
    // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
    virtual bool key_pressed(unsigned int key, int modifiers);

    // This function is called when a keyboard key is down
    // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
    virtual bool key_down(int key, int modifiers);

    // This function is called when a keyboard key is release
    // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
    virtual bool key_up(int key, int modifiers);

    // This function is called when a keyboard key is held down until it repeats
    // - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
    virtual bool key_repeat(int key, int modifiers);

protected:

    // Pointer to the main Viewer class
    Viewer* mViewer;

    // Plugin name
    std::string mName;
};

#endif // VIEWER_PLUGIN_H
