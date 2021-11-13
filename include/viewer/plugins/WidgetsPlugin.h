#ifndef WIDGETS_PLUGIN_H
#define WIDGETS_PLUGIN_H

#include "viewer/ViewerWidget.h"


class Viewer;


class WidgetsPlugin : public ViewerPlugin
{

public:

    WidgetsPlugin();
    ~WidgetsPlugin() override;

    void init(Viewer* viewer) override;
    void shutdown() override;

    // Callbacks
    bool post_load() override;
    bool pre_draw(bool first) override;
    void draw(bool first);
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

    void add(ViewerWidget* widget);

    virtual void reload_font();
    virtual void reload_font(int font_size);

    [[nodiscard]] float menu_scaling() const;
    [[nodiscard]] float pixel_ratio() const;
    [[nodiscard]] float hidpi_scaling() const;
    static float new_pixel_ratio();
    static float new_hidpi_scaling();

private:

    std::vector<ViewerWidget*> widgets;

    // ImGui Context
    ImGuiContext* context_ = nullptr;

    // Ratio between the framebuffer size and the window size.
    // May be different from the hipdi scaling!
    float pixel_ratio_;

    // Hidpi scaling to be used for text rendering.
    float hidpi_scaling_;
};

#endif // WIDGETS_PLUGIN_H
