#ifndef LOADER_PLUGIN_H
#define LOADER_PLUGIN_H

#include "viewer/ViewerWidget.h"
#include "Mesh.h"


class Viewer;


class LoaderPlugin : public ViewerPlugin
{

public:

    LoaderPlugin();
    ~LoaderPlugin() override;

    void init(Viewer* viewer) override;
    void shutdown() override;

    // Callbacks
    bool load(const std::string& filename) override;
    bool unload() override;
    bool save(const std::string& filename) override;
    bool post_load() override;
    bool pre_draw(bool first) override;
    bool post_draw(bool first) override;
    bool post_resize(int w, int h) override;

    // Mouse IO
    bool mouse_down(int button, int modifier) override;
    bool mouse_up(int button, int modifier) override;
    bool mouse_move(int mouse_x, int mouse_y) override;
    bool mouse_scroll(float delta_y) override;

    // Keyboard IO
    bool key_pressed(unsigned int key, int modifiers) override;
    bool key_down(int key, int modifiers) override;
    bool key_up(int key, int modifiers) override;
    bool key_repeat(int key, int modifiers) override;

    // Methods

    [[nodiscard]] Mesh* mesh(int index = -1) const;

private:

    std::vector<Mesh*> mesh_list;
};

#endif // WIDGETS_PLUGIN_H
