#include "viewer/plugins/LoaderPlugin.h"

#include <iostream>


LoaderPlugin::LoaderPlugin() : ViewerPlugin()
{

}

LoaderPlugin::~LoaderPlugin()
{
    for (Mesh* mesh : mesh_list)
    {
        delete mesh;
    }
    mesh_list.clear();
}

void LoaderPlugin::init(Viewer* viewer)
{
    ViewerPlugin::init(viewer);
}

void LoaderPlugin::shutdown()
{

}

bool LoaderPlugin::load(const std::string& filename)
{
    Mesh* mesh = new Mesh();
    if (mesh->load(filename))
    {
        mViewer->append_data(true);
        mesh_list.push_back(mesh);
        assert(mesh_list.size() == mViewer->data_list.size());
        mesh->updateViewerData(mViewer->data_list.back());

        return true;
    }

    delete mesh;
    return false;
}

bool LoaderPlugin::unload()
{
    if (mViewer->selected_data_index < 0 ||
        mViewer->selected_data_index >= mViewer->data_list.size())
    {
        return true;
    }

    delete mesh_list[mViewer->selected_data_index];
    mesh_list.erase(mesh_list.begin() + mViewer->selected_data_index);

    mViewer->data(mViewer->selected_data_index).clear();
    mViewer->data(mViewer->selected_data_index).clear_edges();
    mViewer->erase_data(mViewer->selected_data_index);

    mViewer->selected_data_index = -1;
    return true;
}

bool LoaderPlugin::save(const std::string& filename)
{
    if (mViewer->selected_data_index < 0 ||
        mViewer->selected_data_index >= mViewer->data_list.size())
    {
        return true;
    }

    if (mesh_list[mViewer->selected_data_index])
    {
        return mesh_list[mViewer->selected_data_index]->save(filename);
    }

    return false;
}

bool LoaderPlugin::post_load()
{
    return false;
}

bool LoaderPlugin::pre_draw(bool first)
{
    return false;
}

bool LoaderPlugin::post_draw(bool first)
{
    return false;
}

bool LoaderPlugin::post_resize(int w, int h)
{
    return false;
}

bool LoaderPlugin::mouse_down(int button, int modifier)
{
    return false;
}

bool LoaderPlugin::mouse_up(int button, int modifier)
{
    return false;
}

bool LoaderPlugin::mouse_move(int mouse_x, int mouse_y)
{
    return false;
}

bool LoaderPlugin::mouse_scroll(float delta_y)
{
    return false;
}

bool LoaderPlugin::key_pressed(unsigned int key, int modifiers)
{
    return false;
}

bool LoaderPlugin::key_down(int key, int modifiers)
{
    return false;
}

bool LoaderPlugin::key_up(int key, int modifiers)
{
    return false;
}

bool LoaderPlugin::key_repeat(int key, int modifiers)
{
    return false;
}

Mesh* LoaderPlugin::mesh(int index) const
{
    if (index < 0 || index >= mesh_list.size())
    {
        return nullptr;
    }
    return mesh_list[index];
}
