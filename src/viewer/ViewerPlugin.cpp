#include "viewer/ViewerPlugin.h"


ViewerPlugin::ViewerPlugin()
{
    mViewer = nullptr;
    mName = "dummy";
}

ViewerPlugin::~ViewerPlugin() = default;

void ViewerPlugin::init(Viewer* viewer)
{
    mViewer = viewer;
}

void ViewerPlugin::shutdown()
{

}

bool ViewerPlugin::load(const std::string& filename)
{
    return false;
}

bool ViewerPlugin::unload()
{
    return false;
}

bool ViewerPlugin::save(const std::string& filename)
{
    return false;
}

bool ViewerPlugin::serialize(std::vector<char>& buffer) const
{
    return false;
}

bool ViewerPlugin::deserialize(const std::vector<char>& buffer)
{
    return false;
}

bool ViewerPlugin::post_load()
{
    return false;
}

bool ViewerPlugin::pre_draw(bool first)
{
    return false;
}

bool ViewerPlugin::post_draw(bool first)
{
    return false;
}

bool ViewerPlugin::post_resize(int w, int h)
{
    return false;
}

bool ViewerPlugin::mouse_down(int button, int modifier)
{
    return false;
}

bool ViewerPlugin::mouse_up(int button, int modifier)
{
    return false;
}

bool ViewerPlugin::mouse_move(int mouse_x, int mouse_y)
{
    return false;
}

bool ViewerPlugin::mouse_scroll(float delta_y)
{
    return false;
}

bool ViewerPlugin::key_pressed(unsigned int key, int modifiers)
{
    return false;
}

bool ViewerPlugin::key_down(int key, int modifiers)
{
    return false;
}

bool ViewerPlugin::key_up(int key, int modifiers)
{
    return false;
}

bool ViewerPlugin::key_repeat(int key, int modifiers)
{
    return false;
}
